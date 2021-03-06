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

package Bio::EnsEMBL::Tark::HGNC;

use Moose;
with 'MooseX::Log::Log4perl';

use Bio::EnsEMBL::Tark::DB;
use Bio::EnsEMBL::Tark::FileHandle;

use Try::Tiny;

has 'query' => (
  traits  => ['Hash'],
  is      => 'rw',
  isa     => 'HashRef',
  default => sub { {} },
  handles => {
    set_query     => 'set',
    get_query     => 'get',
    delete_query  => 'delete',
    clear_queries => 'clear',
    fetch_keys    => 'keys',
    fetch_values  =>'values',
    query_pairs   => 'kv',
  },
);

has session => (
  is  => 'rw',
  isa => 'Bio::EnsEMBL::Tark::DB',
);


=head2 BUILD
  Description: Initialise the creation of the object
  Returntype : undef
  Exceptions : none
  Caller     : general

=cut

sub BUILD {
  my ($self) = @_;

  $self->log()->info('Initializing HGNC loader');

  # Attempt a connection to the database
  my $dbh = $self->session->dbh();

  # Setup the insert queries
  my $insert_gene_name_sql = (<<'SQL');
    INSERT INTO gene_names (external_id, name, source, primary_id, session_id)
    VALUES (?, ?, 'HGNC', ?, ?)
SQL

  my $sth = $dbh->prepare( $insert_gene_name_sql );
  if ( $sth->err ) {
    $self->log->logdie('Error creating gene name insert: ' . $sth->errstr);
  }
  $self->set_query('hgnc' => $sth);

  my $get_ids_sql = (<<'SQL');
    SELECT gene_id, assembly_id FROM gene WHERE stable_id = ?
SQL

  $sth = $dbh->prepare( $get_ids_sql );
  if ( $sth->err ) {
    $self->log->logdie('Error creating gene select: ' . $sth->errstr);
  }
  $self->set_query('gene' => $sth);

  return;
} ## end sub BUILD


=head2 flush_hgnc
  Description: Truncate the gene_names table
  Returntype : undef
  Exceptions : none
  Caller     : general

=cut

sub flush_hgnc {
  my $self = shift;

  my $dbh = $self->session->dbh();

  $self->log()->info('Truncating gene names table');
  $dbh->do('TRUNCATE gene_names');

  return;
} ## end sub flush_hgnc


=head2 load_hgnc
  Arg [1]    : $hgnc_file : string
  Description: Extract HGNC names and aliases from the HGNC dump files
  Returntype : undef
  Exceptions : none
  Caller     : general
  Notes      : Col 1: hgnc_id
               Col 2: symbol/name
               Col 9: alias_symbols
               Col 20: ensembl_gene_id
              (counting from 1)

=cut

sub load_hgnc {
  my $self = shift;
  my $hgnc_file = shift;

  $self->log()->info('Starting HGNC load');

  my $in_fh;
  if($hgnc_file) {
    $self->log()->info("Using HGNC file $hgnc_file");

    my $file_handle = Bio::EnsEMBL::Tark::FileHandle->new();
    $in_fh = $file_handle->get_file_handle($hgnc_file);
  } else {
    $in_fh = *STDIN;
  }

  my $get_gene        = $self->get_query('gene');
  my $insert_hgnc     = $self->get_query('hgnc');

  while(<$in_fh>) {

    chomp;

    my @hgnc_line = split '\t';
    if ( $hgnc_line[0] !~ /^HGNC:/ ) {
      next;
    }

    # Insert the hgnc symbol
    $insert_hgnc->execute($hgnc_line[0], $hgnc_line[1], 1, $self->session->session_id);

    # Add any synomyms
    if ( !$hgnc_line[8] ) {
      next;
    }

    $hgnc_line[8] =~ s/^"//;
    $hgnc_line[8] =~ s/"$//;

    my @aliases = split '\|', $hgnc_line[8];
    foreach my $alias (@aliases) {
      try {
        $insert_hgnc->execute($hgnc_line[0], $alias, 0, $self->session->session_id);
      } catch {
        $self->log->warn("Failed to load alias $alias for $hgnc_line[0]");
      }
    }
  }

  return;
} ## end sub load_hgnc

1;
