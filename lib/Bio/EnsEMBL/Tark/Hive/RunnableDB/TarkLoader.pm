=head1 LICENSE

See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

=cut

package Bio::EnsEMBL::Tark::Hive::RunnableDB::TarkLoader;

use strict;
use warnings;
use Carp;

use base ('Bio::EnsEMBL::Hive::Process');


sub param_defaults {
  return {
    'column_names'      => 0,
    'delimiter'         => undef,
    'randomize'         => 0,
    'step'              => 0,
    'contiguous'        => 0,
    'key_column'        => 0,

    # this parameter is no longer supported and should stay at 0
    'input_id'          => 0,

    'inputlist'         => undef,
    'inputfile'         => undef,
    'inputquery'        => undef,
    'inputcmd'          => undef,

    'fan_branch_code'   => 2,

    # Boolean. When true, the command will be run with "bash -o pipefail -c $cmd".
    # Useful to capture errors in a command that contains pipes
    'use_bash_pipefail' => 0,
  };
}



=head2 run
  Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process
                that is used to perform the main bulk of the job (minus input and
                output).
  param
    column_names : Controls the column names that come out of the parser:
                     * 0 = "no names",
                     * 1 = "parse names from data",
                     * arrayref = "take names from this array"
    delimiter    : If you set it your lines in file/cmd mode will be split into
                   columns that you can use individually when constructing the
                   input_id_template hash.
    randomize    : Shuffles the rows before creating jobs - can sometimes lead to
                   better overall performance of the pipeline. Doesn't make any
                   sence for minibatches (step>1).
    step         : The requested size of the minibatch (1 by default). The real
                   size of a range may be smaller than the requested size.
    contiguous   : Whether the key_column range of each minibatch should be
                   contiguous (0 by default).
    key_column   : If every line of your input is a list (it happens, for example,
                   when your SQL returns multiple columns or you have set the
                   'delimiter' in file/cmd mode) this is the way to say which
                   column is undergoing 'ranging'

  The following 4 parameters are mutually exclusive and define the source of ids
  for the jobs:
    inputlist    : The list is explicitly given in the parameters, can be
                   abbreviated:
                     * 'inputlist' => ['a'..'z']
    inputfile    : The list is contained in a file whose name is supplied as
                   parameter:
                     * 'inputfile' => 'myfile.txt'
    inputquery   : The list is generated by an SQL query (against the production
                   database by default):
                     * 'inputquery' => 'SELECT object_id FROM object WHERE x=y'
    inputcmd     : The list is generated by running a system command:
                     * 'inputcmd' => 'find /tmp/big_directory -type f'
=cut

sub run {
  my $self = shift @_;

  # can be 0 (no names), 1 (names from data) or an arrayref (names from this array)
  my $column_names = $self->param('column_names');
  my $delimiter    = $self->param('delimiter');

  my $randomize    = $self->param('randomize');

      # minibatching-related:
  my $step         = $self->param('step');
  my $contiguous   = $self->param('contiguous');
  my $key_column   = $self->param('key_column');

  my $inputlist    = $self->param('inputlist');
  my $inputfile    = $self->param('inputfile');
  my $inputquery   = $self->param('inputquery');
  my $inputcmd     = $self->param('inputcmd');

  my $parse_column_names = $column_names && (ref($column_names) ne 'ARRAY');

  my ($rows, $column_names_from_data) =
      $inputlist    ? $self->_get_rows_from_list(  $inputlist  )
    : $inputquery   ? $self->_get_rows_from_query( $inputquery )
    : $inputfile    ? $self->_get_rows_from_open(  $inputfile  , '<', $delimiter, $parse_column_names )
    : $inputcmd     ? $self->_get_rows_from_open( ($self->param('use_bash_pipefail') ? 'set -o pipefail; ': '').$inputcmd, '-|', $delimiter, $parse_column_names )
    : confess "range of values should be defined by setting 'inputlist', 'inputquery', 'inputfile' or 'inputcmd'";

  if(
    $column_names_from_data                         # column data is available
  and ( defined($column_names) ? (ref($column_names) ne 'ARRAY') : 1 )    # and is badly needed
  ) {
    $column_names = $column_names_from_data;
  }
  # after this point $column_names should either contain a list or be false

  if( $self->param('input_id') ) {
    confess "'input_id' is no longer supported, please reconfigure as the input_id_template of the dataflow_rule";
  }

  if($randomize) {
      _fisher_yates_shuffle_in_place($rows);
  }

  my $output_ids = $step
    ? $self->_substitute_minibatched_rows($rows, $column_names, $step, $contiguous, $key_column)
    : $self->_substitute_rows($rows, $column_names);

  $self->param('output_ids', $output_ids);
}


=head2 write_output
  Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                Here we rely on the dataflow mechanism to create jobs.
  param('fan_branch_code'): defines the branch where the fan of jobs is created (2 by default).
=cut

sub write_output {  # nothing to write out, but some dataflow to perform:
    my $self = shift @_;

    my $output_ids              = $self->param('output_ids');
    my $fan_branch_code         = $self->param('fan_branch_code');

        # "fan out" into fan_branch_code:
    $self->dataflow_output_id($output_ids, $fan_branch_code);
}


################################### main functionality starts here ###################




1;