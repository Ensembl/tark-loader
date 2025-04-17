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

=head1 DESCRIPTION

  A pipeline for the loading of the Tark DB from a Core db.

=head1 SYNOPSIS

  # Prepare the pipeline
  init_pipeline.pl Bio::EnsEMBL::Tark::Hive::PipeConfig::TarkLoader_conf
    --tark_host $TARK_HOST
    --tark_port $TARK_PORT
    --tark_user $TARK_USER
    --tark_pass $TARK_PASS
    --tark_db $TARK_SPP_DB
    --core_host $CORE_HOST
    --core_port $CORE_PORT
    --core_user $CORE_USER
    --core_pass $CORE_PASS
    --core_dbname $ENS_SPP_CORE_DB
    --host $HIVE_HOST
    --port $HIVE_PORT
    --user $HIVE_USER
    --password $HIVE_PASS
    --pipeline_name test_hive_1234
    --species homo_sapiens
    --tag_block release
    --tag_shortname 84
    --tag_release_date 2014-10-01
    --tag_description 'Ensembl release 84'
    --tag_feature_type all
    --tag_version 1
    --block_size 1000

=cut

package Bio::EnsEMBL::Tark::Hive::PipeConfig::TarkLoader_conf;

use strict;
use warnings;

# All Hive databases configuration files should inherit from HiveGeneric, directly
# or indirectly
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

# Allow this particular config to use conditional dataflow and INPUT_PLUS
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

use Bio::EnsEMBL::Tark::Hive::PipeConfig::SQL;

use Bio::EnsEMBL::DBSQL::DBAdaptor;



=head2 default_options
    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines four options:
                    o('concurrent_jobs')   defines how many tables can be worked on in parallel
=cut

sub default_options {
    my ($self) = @_;

    return {
      %{$self->SUPER::default_options(@_)},

      'core_db'  => {
        -host   => $self->o( 'core_host' ),
        -port   => $self->o( 'core_port' ),
        -user   => $self->o( 'core_user' ),
        -pass   => $self->o( 'core_pass' ),
        -dbname => $self->o( 'core_dbname' ),
        -driver => 'mysql',
      },

      exclude_source         => q{},
      include_source         => q{},
      tag_previous_shortname => q{},
      naming_consortium      => q{},
      add_consortium_prefix  => 0,
    };
} ## end sub default_options

=head2 pipeline_analyses
  Description : Implements pipeline_analyses() interface method of
                Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the
                structure of the pipeline: analyses, jobs, rules, etc.
                Here it defines two analyses:
                  * 'job_factory' which uses the JobFactory runnable to optain a
                    chunked set of lists of gene_ids from the Core db.
                  * 'TarkLoader' which takes a chunk produced by the JobFactory
                    and loads the listed genes into the Tark db.

=cut

sub resource_classes {
    my ($self) = @_;
    return {
        'default' => { SLURM => '--partition=production --time=1-00:00:00 --mem=4000' },
        'default_slurm' => { SLURM => '--partition=production --time=1-00:00:00 --mem=4000' },
        'default_16GB'    => { SLURM => '--partition=production --time=1-00:00:00 --mem=16000' },
    };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'generate_sql',
      -module     => 'Bio::EnsEMBL::Tark::Hive::RunnableDB::GeneSetSQL',
      -parameters => {
        block_size => $self->o('block_size'),
        exclude_source => $self->o('exclude_source'),
        include_source => $self->o('include_source'),
      },
      -input_ids => [
        {
          'target_db'  => $self->dbconn_2_url( 'core_db' ),
        }
      ],
      -flow_into  => {
        '2->A' => { 'generate_sql_params' => INPUT_PLUS() },
        'A->1' => { 'tark_report'         => INPUT_PLUS() }
      },
    },

    {
      -logic_name => 'generate_sql_params',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        'db_conn'      => $self->dbconn_2_url( 'core_db' ),
        'column_names' => [ 'gene_id_list' ],
        'inputquery'   => '#sql#',
        },
      -flow_into  => { 2 => { 'load_tark' => INPUT_PLUS() } },
    },

    {
      -logic_name => 'load_tark',
      -module     => 'Bio::EnsEMBL::Tark::Hive::RunnableDB::TarkLoader',
      -parameters => {

        # Worker block size params
        tag_block        => $self->o( 'tag_block' ),
        tag_shortname    => $self->o( 'tag_shortname' ),
        tag_description  => $self->o( 'tag_description' ),
        tag_feature_type => $self->o( 'tag_feature_type' ),
        tag_version      => $self->o( 'tag_version' ),
        tag_release_date => $self->o( 'tag_release_date' ),

        # Species name
        species     => $self->o( 'species' ),

        # Core db params
        host => $self->o( 'core_host' ),
        port => $self->o( 'core_port' ),
        user => $self->o( 'core_user' ),
        pass => $self->o( 'core_pass' ),
        db   => $self->o( 'core_dbname' ),

        # Tark db params
        tark_host => $self->o( 'tark_host' ),
        tark_port => $self->o( 'tark_port' ),
        tark_user => $self->o( 'tark_user' ),
        tark_pass => $self->o( 'tark_pass' ),
        tark_db   => $self->o( 'tark_db' ),

        source_name           => $self->o( 'source_name' ),
        naming_consortium     => $self->o( 'naming_consortium' ),
        add_consortium_prefix => $self->o( 'add_consortium_prefix' ),
      },
      -analysis_capacity => 50,
    },

    {
      -logic_name => 'tark_report',
      -module     => 'Bio::EnsEMBL::Tark::Hive::RunnableDB::TarkLoaderReport',
      -parameters => {

        # Core db params
        host => $self->o( 'core_host' ),
        port => $self->o( 'core_port' ),
        user => $self->o( 'core_user' ),
        pass => $self->o( 'core_pass' ),
        db   => $self->o( 'core_dbname' ),

        # Tark db params
        tark_host => $self->o( 'tark_host' ),
        tark_port => $self->o( 'tark_port' ),
        tark_user => $self->o( 'tark_user' ),
        tark_pass => $self->o( 'tark_pass' ),
        tark_db   => $self->o( 'tark_db' ),

        # Include/Exclude db entries by source
        exclude_source => $self->o( 'exclude_source' ),
        include_source => $self->o( 'include_source' ),

        # Worker parameters
        source_name            => $self->o( 'source_name' ),
        tag_shortname          => $self->o( 'tag_shortname' ),
        tag_release_date => $self->o( 'tag_release_date' ),
        tag_previous_shortname => $self->o( 'tag_previous_shortname' ),
        report                 => $self->o( 'report' ),
      }
    }
  ];
} ## end sub pipeline_analyses

1;
