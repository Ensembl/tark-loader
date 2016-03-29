#!/usr/bin/env perl


=head1 LICENSE
Copyright 2016 EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Tark::SpeciesLoader;
use Bio::EnsEMBL::Registry;

my $dbuser; my $dbpass; my $dbhost; my $database; my $dbport = 3306;
my $species;
my $ensdbhost = 'mysql-ensembl-mirror.ebi.ac.uk';
my $ensdbport = 4240;

Log::Log4perl->easy_init($DEBUG);

get_options();

my $loader = Bio::EnsEMBL::Tark::SpeciesLoader->new( dsn => "DBI:mysql:database=$database;host=$dbhost;port=$dbport",
						     dbuser => $dbuser,
						     dbpass => $dbpass );

# Connect to the Ensembl Registry to access the databases
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $ensdbhost,
    -port => $ensdbport,
    -user => 'anonymous',
    -db_version => '84'
    );

my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, "core" );

my $session_id = $loader->start_session("Test client");

#eval {
    $loader->load_species($dba, $session_id);
#};
#if($@) {
#    $loader->abort_session($session_id);
#    die "Error loading species: $@";
#}

$loader->end_session($session_id);


sub get_options {
    my $help;

    GetOptions(
	"dbuser=s"               => \$dbuser,
	"dbpass=s"               => \$dbpass,
	"dbhost=s"               => \$dbhost,
	"database=s"             => \$database,
	"dbport=s"               => \$dbport,
        "species=s"              => \$species,
	"enshost=s"              => \$ensdbhost,
	"ensport=s"              => \$ensdbport,
        "help"                   => \$help,
        );
    
    if ($help) {
        exec('perldoc', $0);
    }

}
