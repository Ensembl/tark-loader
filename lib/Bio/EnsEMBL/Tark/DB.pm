=head1 LICENSE

Copyright 2015 EMBL-European Bioinformatics Institute

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

use warnings;
use strict;
use DBI;
use Carp;
use Digest::SHA1 qw(sha1);

package Bio::EnsEMBL::Tark::DB;

use MooseX::Singleton;
with 'MooseX::Log::Log4perl';

my $singleton;

has 'dsn' => ( is => 'ro', isa => 'Str' );

has 'dbuser' => ( is => 'ro', isa => 'Str' );

has 'dbpass' => ( is => 'ro', isa => 'Str' );

has 'session_id' => ( is => 'rw', isa => 'Int', default => 0 );

sub dbh {
    my $self = shift;

    return DBI->connect_cached( $self->dsn, $self->dbuser, $self->dbpass )
	or $self->log()->die("Error connecting to " . $self->dsn . ": ". $DBI::errstr);

}

# Join an array of values with a ':' delimeter and find a sha1 checksum of it

sub checksum_array {
    my ($self, @values) = @_;

    return Digest::SHA1::sha1( join(':', grep { defined } @values) );
}

sub start_session {
    my $self = shift;
    my $client_name = shift;

    my $dbh = $self->dbh();
    my $sth = $dbh->prepare("INSERT INTO session (client_id, status) VALUES(?, 1)");
    $sth->execute($client_name) or
		$self->log->logdie("Error inserting session: $DBI::errstr");

    $self->session_id($sth->{mysql_insertid});

    $self->log->info("Starting session " . $self->session_id);

#    $dbh->do("SET FOREIGN_KEY_CHECKS = 0");
#    $dbh->do("SET UNIQUE_CHECKS = 0");
#    $dbh->do("SET SESSION tx_isolation='READ-UNCOMMITTED'");
#    $dbh->do("SET sql_log_bin = 0");

    return $self->session_id;
}

sub end_session {
    my $self = shift;

    return unless $self->session_id;

    my $dbh = $self->dbh();
    $dbh->do("UPDATE session SET status = 2 WHERE session_id = ?", undef, $self->session_id);

    $self->session_id(0);

#    $dbh->do("SET UNIQUE_CHECKS = 1");
#    $dbh->do("SET FOREIGN_KEY_CHECKS = 1");
#    $dbh->do("SET SESSION tx_isolation='READ-REPEATABLE'");

}

sub abort_session {
    my $self = shift;

    return unless $self->session_id;

    my $dbh = $self->dbh();
    $dbh->do("UPDATE session SET status = 3 WHERE session_id = ?", undef, $self->session_id);

    $self->session_id(0);

    $dbh->do("SET UNIQUE_CHECKS = 1");
    $dbh->do("SET FOREIGN_KEY_CHECKS = 1");
}

1;
