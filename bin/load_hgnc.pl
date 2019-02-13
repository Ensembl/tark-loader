#!/usr/bin/env perl


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

$|++;

use warnings;
use strict;

use Log::Log4perl qw(:easy);
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Tark::HGNC;
use Bio::EnsEMBL::Tark::DB;

my ( $dbuser, $dbpass, $dbhost, $database, $hgnc_file, $flush_names );
my $dbport = 3306;

Log::Log4perl->easy_init( $DEBUG );

get_options();

my $db = Bio::EnsEMBL::Tark::DB->new(
  config => {
    driver => 'mysql',
    host   => $dbhost,
    port   => $dbport,
    user   => $dbuser,
    pass   => $dbpass,
    db     => $database,
  }
);

my $session_id = $db->start_session( 'HGNC loader' );

my $loader = Bio::EnsEMBL::Tark::HGNC->new(
  session => $db
);

if($flush_names) {
  $loader->flush_hgnc();
}

$loader->load_hgnc( $hgnc_file );

$db->end_session( $session_id );


=head2 get_options
  Description:
  Returntype :
  Exceptions : none
  Caller     : general

=cut

sub get_options {
  my $help;

  GetOptions(
    # 'config=s' => \$config_file,
    'dbuser=s'   => \$dbuser,
    'dbpass=s'   => \$dbpass,
    'dbhost=s'   => \$dbhost,
    'database=s' => \$database,
    'dbport=s'   => \$dbport,
    'hgnc=s'     => \$hgnc_file,
    'flush'      => \$flush_names,
    'help'       => \$help,
  );

  if ($help) {
    exec 'perldoc', $0;
  }
} ## end sub get_options
