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

package Bio::EnsEMBL::Tark::SpeciesLoader;

use Moose;
with 'MooseX::Log::Log4perl';

use DBI;

use Bio::EnsEMBL::Tark::DB;
use Bio::EnsEMBL::Tark::Tag;
use Bio::EnsEMBL::Tark::Utils;


has 'dsn' => ( is => 'ro', isa => 'Str' );

has 'dbuser' => ( is => 'ro', isa => 'Str' );

has 'dbpass' => ( is => 'ro', isa => 'Str' );

has 'query' => (
  traits  => ['Hash'],
  is      => 'rw',
  isa     => 'HashRef',
  default => sub { {} },
  handles => {
    set_query     => 'set',
    get_insert    => 'get',
    delete_query  => 'delete',
    clear_queries => 'clear',
    fetch_keys    => 'keys',
    fetch_values  => 'values',
    query_pairs   => 'kv',
  },
);

has session => (
  is  => 'rw',
  isa => 'Bio::EnsEMBL::Tark::DB',
);

has tag_config => (
  is  => 'rw',
  isa => 'Bio::EnsEMBL::Tark::TagConfig'
);

has gene_id_list => ( is => 'ro', isa => 'ArrayRef' );

has naming_consortium => ( is => 'ro', isa => 'Str' );
has add_name_prefix => ( is => 'ro', isa => 'Int' );


=head2 BUILD
  Description: Initialise the creation of the prepared statements
  Returntype : undef
  Exceptions : none
  Caller     : general

=cut

sub BUILD {
  my ( $self, $args ) = @_;

  $self->log()->info('Initializing species loader');

  # Attempt a connection to the database
  my $dbh = $self->session->dbh();

  # Setup the insert queries

  # INSERT genome
  my $genome_sql = (<<'SQL');
    INSERT INTO genome (name, tax_id, session_id)
    VALUES (?, ?, ?)
    ON DUPLICATE KEY UPDATE genome_id=LAST_INSERT_ID(genome_id)
SQL

  my $sth = $dbh->prepare( $genome_sql ) or
    $self->log->logdie("Error creating genome insert: $DBI::errstr");
  $self->set_query('genome' => $sth);


  # INSERT assembly
  my $assembly_sql = (<<'SQL');
    INSERT INTO assembly (genome_id, assembly_name, session_id)
    VALUES (?, ?, ?)
    ON DUPLICATE KEY UPDATE assembly_id=LAST_INSERT_ID(assembly_id)
SQL

  $sth = $dbh->prepare( $assembly_sql );
  $self->set_query('assembly' => $sth);


  # INSERT assembly_alias
  my $assembly_alias_sql = (<<'SQL');
    INSERT INTO assembly_alias (genome_id, assembly_id, alias, session_id)
    VALUES (?, ?, ?, ?)
    ON DUPLICATE KEY UPDATE assembly_id=LAST_INSERT_ID(assembly_id)
SQL

  $sth = $dbh->prepare( $assembly_alias_sql );
  $self->set_query('assembly_alias' => $sth);

  
  # INSERT gene
  # populate biotype for gene
  my $gene_sql = (<<'SQL');
    INSERT INTO gene (
      stable_id, stable_id_version, assembly_id, loc_region, loc_start, loc_end,
      loc_strand, loc_checksum, name_id, gene_checksum, session_id, biotype)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON DUPLICATE KEY UPDATE gene_id=LAST_INSERT_ID(gene_id)
SQL
  $sth = $dbh->prepare( $gene_sql );
  $self->set_query('gene' => $sth);


  # INSERT transcript_gene
  my $transcript_gene_sql = (<<'SQL');
    INSERT INTO transcript_gene (gene_id, transcript_id, session_id)
    VALUES (?, ?, ?)
    ON DUPLICATE KEY UPDATE gene_transcript_id=LAST_INSERT_ID(gene_transcript_id)
SQL
  $sth = $dbh->prepare( $transcript_gene_sql );
  $self->set_query('transcript_gene' => $sth);


  # INSERT transcript
  # populate biotype for transcript
  my $transcript_sql = (<<'SQL');
    INSERT INTO transcript (
      stable_id, stable_id_version, assembly_id, loc_region, loc_start, loc_end,
      loc_strand, loc_checksum, transcript_checksum, exon_set_checksum,
      seq_checksum, session_id, biotype, three_prime_utr_start, three_prime_utr_end,
      three_prime_utr_seq, three_prime_utr_checksum, five_prime_utr_start,
      five_prime_utr_end, five_prime_utr_seq, five_prime_utr_checksum)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    ON DUPLICATE KEY UPDATE transcript_id=LAST_INSERT_ID(transcript_id)
SQL
  $sth = $dbh->prepare( $transcript_sql );
  $self->set_query('transcript' => $sth);


  # INSERT exon_transcript
  my $exon_transcript_sql = (<<'SQL');
    INSERT INTO exon_transcript (transcript_id, exon_id, exon_order, session_id)
    VALUES (?, ?, ?, ?)
    ON DUPLICATE KEY UPDATE exon_transcript_id=LAST_INSERT_ID(exon_transcript_id)
SQL
  $sth = $dbh->prepare( $exon_transcript_sql );
  $self->set_query('exon_transcript' => $sth);


  # INSERT exon
  my $exon_sql = (<<'SQL');
    INSERT INTO exon (
      stable_id, stable_id_version, assembly_id, loc_region, loc_start, loc_end,
      loc_strand, loc_checksum, exon_checksum, seq_checksum, session_id)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON DUPLICATE KEY UPDATE exon_id=LAST_INSERT_ID(exon_id)
SQL
  $sth = $dbh->prepare( $exon_sql );
  $self->set_query('exon' => $sth);


  # INSERT translation
  my $translation_sql = (<<'SQL');
    INSERT INTO translation (
      stable_id, stable_id_version, assembly_id, loc_region, loc_start, loc_end,
      loc_strand, loc_checksum, translation_checksum, seq_checksum, session_id)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON DUPLICATE KEY UPDATE translation_id=LAST_INSERT_ID(translation_id)
SQL
  $sth = $dbh->prepare( $translation_sql );
  $self->set_query('translation' => $sth);


  # INSERT translation_transcript
  my $translation_transcript_sql = (<<'SQL');
    INSERT INTO translation_transcript (transcript_id, translation_id, session_id)
    VALUES (?, ?, ?)
    ON DUPLICATE KEY UPDATE transcript_translation_id=LAST_INSERT_ID(transcript_translation_id)
SQL
  $sth = $dbh->prepare( $translation_transcript_sql );
  $self->set_query('translation_transcript' => $sth);


  # INSERT sequence
  my $sequence_sql = (<<'SQL');
    INSERT IGNORE INTO sequence (seq_checksum, sequence, session_id)
    VALUES (?, ?, ?)
SQL
  $sth = $dbh->prepare( $sequence_sql );
  $self->set_query('sequence' => $sth);

  return;
} ## end sub BUILD


=head2 load_species
  Description: Load a species worth of genomic featues into the Tark DB from an
               ensmbl core db.
  Returntype : undef
  Exceptions : none
  Caller     : general

=cut

sub load_species {
  my ( $self, $dba, $source_name ) = @_;

  my $session_id = $self->session->session_id;

  my $naming_consortium = $self->naming_consortium;
  my $add_name_prefix   = $self->add_name_prefix;

  $source_name = defined($source_name) ? $source_name : 'Ensembl';

  $self->log->info( 'Starting loading process' );

  my $mc = $dba->get_MetaContainer();
  my $species = $mc->get_production_name();
  $self->log->info( "Storing genome for $species" );

  my $sth = $self->get_insert( 'genome' );
  $sth->execute( $species, $mc->get_taxonomy_id() + 0, $session_id ) or
    $self->log->logdie( "Error inserting genome: $DBI::errstr" );

  my $genome_id = $sth->{ mysql_insertid };
  $sth = $self->get_insert( 'assembly' );

  my $assembly_accession = $mc->single_value_by_key( 'assembly.accession' );
  my $assembly_name = $mc->single_value_by_key( 'assembly.default' );

  $self->log->info( "Inserting assembly $assembly_name" );
  $sth->execute( $genome_id, $assembly_name, $session_id ) or
    $self->log->logdie( "Error inserting assembly: $DBI::errstr" );

  my $assembly_id = $sth->{ mysql_insertid };
  $self->log->info( "Using assembly_id $assembly_name (id: $assembly_id)" );

  $sth = $self->get_insert( 'assembly_alias' );
  $sth->execute( $genome_id, $assembly_id, $assembly_accession, $session_id ) or
    $self->log->logdie( "Error inserting assembly_alias: $DBI::errstr" );


  # Initialize the tags we'll be using
  my $tag = Bio::EnsEMBL::Tark::Tag->new(
    session => $self->session,
    config  => $self->tag_config
  );
  $tag->init_tags( $assembly_id );

  my $session_pkg = {
    session_id  => $session_id,
    genome_id   => $genome_id,
    assembly_id => $assembly_id,
    assembly_name => $assembly_name
  };

  my $iter;
  if ( defined $self->gene_id_list ) {
    # Fetch and load genes within a defined ID range
    my $ga   = $dba->get_GeneAdaptor();
    foreach my $current_gene ( @{ $self->gene_id_list } ) {
      my $gene = $ga->fetch_by_dbID( $current_gene );
      if ( $gene ) {
        $self->log->debug( 'Loading gene ' . $gene->{stable_id} );
        $self->_load_gene(
          $gene, $session_pkg, $source_name, $tag, $naming_consortium,
          $add_name_prefix
        );
      }
    }
  } else {
    # Fetch a gene iterator and cycle through loading the genes
    $iter = $self->genes_to_metadata_iterator( $dba, $source_name );

    while ( my $gene = $iter->next() ) {
      $self->log->debug( 'Loading gene ' . $gene->{stable_id} );
      $self->_load_gene(
        $gene, $session_pkg, $source_name, $tag, $naming_consortium,
        $add_name_prefix
      );
    }
  }

  $self->log->info( 'Completed dumping genes for ' . $species );

  $self->log->info( 'Tagging sets for ' . $species );
  $tag->checksum_sets();

  $self->log->info( 'Completed tagging sets for ' . $species );

  return;
} ## end sub load_species


=head2 _load_gene
  Arg [1]    : $gene        : Bio::EnsEMBL::Gene
  Arg [2]    : $session_pkg : Bio::EnsEMBL::Tark::DB
  Arg [4]    : $source_name : string
  Arg [4]    : $tag         : Bio::EnsEMBL::Tark::Tag
  Description: For each gene generate the matching checksum and initiate the
               creation of checksums for the exons, transcripts and translations
  Returntype : undef
  Exceptions : none
  Caller     : general

=cut

sub _load_gene {
  my ( $self, $gene, $session_pkg, $source_name, $tag, $naming_consortium, $add_name_prefix ) = @_;

  my @loc_pieces = (
    $session_pkg->{assembly_id}, $gene->seq_region_name(),
    $gene->seq_region_start(), $gene->seq_region_end(),
    $gene->seq_region_strand(),
  );

  my @loc_pieces_checksum = (
    $session_pkg->{assembly_name}, $gene->seq_region_name(),
    $gene->seq_region_start(), $gene->seq_region_end(),
    $gene->seq_region_strand(),
  );

  my $utils = Bio::EnsEMBL::Tark::Utils->new();
  my $loc_checksum = $utils->checksum_array( @loc_pieces_checksum );

  my $name_id = undef;
  if ( $naming_consortium ) {
    $name_id = $self->_fetch_name_id($gene, $naming_consortium, $add_name_prefix);
  }

  my $gene_checksum = $utils->checksum_array(
    @loc_pieces, $name_id, $gene->stable_id(), $gene->version(), $gene->biotype()
  );

  my $sth = $self->get_insert('gene');
  # populate biotype for gene
  $sth->execute(
    $gene->stable_id(), $gene->version(), @loc_pieces, $loc_checksum,
    $name_id, $gene_checksum, $session_pkg->{session_id}, $gene->biotype()
  ) or  $self->log->logdie("Error inserting gene: $DBI::errstr");

  my $gene_id = $sth->{mysql_insertid};

  # Apply tags to feature we've just inserted
  $tag->tag_feature($gene_id, 'gene');

  my $exons = {};
  $session_pkg->{gene_id} = $gene_id;
  for my $transcript ( @{ $gene->get_all_Transcripts() } ) {

    my @exon_checksums;
    my @exon_ids;
    for my $exon (@{ $transcript->get_all_Exons() }) {
      my ($exon_id, $exon_checksum) = $self->_load_exon( $exon, $session_pkg, $tag );
      push @exon_checksums, $exon_checksum;
      push @exon_ids, $exon_id;
    }

    if( @exon_checksums ) {
      $session_pkg->{exon_set_checksum} =  $utils->checksum_array( @exon_checksums );
    }

    my $transcript_id = $self->_load_transcript( $transcript, $session_pkg, $tag );

    my $exon_order = 1;
    for my $exon_id (@exon_ids) {
      $sth = $self->get_insert('exon_transcript');
      $sth->execute(
        $transcript_id, $exon_id, $exon_order, $session_pkg->{session_id}
      ) or $self->log->logdie("Error inserting exon_transcript: $DBI::errstr");
      $exon_order++;
    }

    $session_pkg->{transcript_id} = $transcript_id;
    $session_pkg->{transcript} = $transcript;

    my $translation = $transcript->translation();
    if ( defined $translation ) {
      $self->_load_translation( $translation, $session_pkg, $tag );
    }

    delete $session_pkg->{transcript_id};
    delete $session_pkg->{transcript};
  }

  return;
} ## end sub _load_gene


=head2 _load_transcript
  Arg [1]    : $transcript  : Bio::EnsEMBL::Transcript
  Arg [2]    : $session_pkg : Bio::EnsEMBL::Tark::DB
  Arg [3]    : $tag         : Bio::EnsEMBL::Tark::Tag
  Description: For each transcript generate and load the checksum
  Returntype : $transcript_id : integer
  Exceptions : none
  Caller     : general

=cut

sub _load_transcript {
  my ( $self, $transcript, $session_pkg, $tag ) = @_;

  # Insert the sequence and get back the checksum
  my $seq_obj = $transcript->seq();
  my $seq_checksum;
  if($seq_obj) {
    $seq_checksum = $self->_insert_sequence(
      $seq_obj->seq(), $session_pkg->{session_id}
    );
  }

  my @loc_pieces = (
    $session_pkg->{assembly_id}, $transcript->seq_region_name(),
    $transcript->seq_region_start(), $transcript->seq_region_end(),
    $transcript->seq_region_strand(),
  );

  my @loc_pieces_checksum = (
    $session_pkg->{assembly_name}, $transcript->seq_region_name(),
    $transcript->seq_region_start(), $transcript->seq_region_end(),
    $transcript->seq_region_strand(),
  );

  my $utils = Bio::EnsEMBL::Tark::Utils->new();
  my $loc_checksum = $utils->checksum_array( @loc_pieces_checksum );
  my $transcript_checksum = $utils->checksum_array(
    $loc_checksum, $transcript->stable_id(), $transcript->version(), $transcript->biotype(),
    (
      $session_pkg->{exon_set_checksum} ? $session_pkg->{exon_set_checksum} : undef
    ), $seq_checksum
  );

  my $sth = $self->get_insert('transcript');
  # populate biotype for transcript
  $sth->execute(
    $transcript->stable_id(), $transcript->version(), @loc_pieces, $loc_checksum,
    $transcript_checksum, (
      $session_pkg->{exon_set_checksum} ? $session_pkg->{exon_set_checksum} : undef
    ), $seq_checksum, $session_pkg->{session_id}, $transcript->biotype()
  ) or $self->log->logdie("Error inserting transcript: $DBI::errstr");
  my $transcript_id = $sth->{mysql_insertid};

  $sth = $self->get_insert('transcript_gene');
  $sth->execute(
    $session_pkg->{gene_id}, $transcript_id, $session_pkg->{session_id}
  ) or $self->log->logdie("Error inserting transcript_gene: $DBI::errstr");

  # Apply tags to feature we've just inserted
  $tag->tag_feature($transcript_id, 'transcript');

  return $transcript_id;
} ## end sub _load_transcript


=head2 _load_exon
  Arg [1]    : $exon        : Bio::EnsEMBL::Exon
  Arg [2]    : $session_pkg : Bio::EnsEMBL::Tark::DB
  Arg [3]    : $tag         : Bio::EnsEMBL::Tark::Tag
  Description: For an exon generate the matching checksum
  Returntype : Array ($exon_id, $exon_checksum)
  Exceptions : none
  Caller     : general

=cut

sub _load_exon {
  my ( $self, $exon, $session_pkg, $tag ) = @_;

  # Insert the sequence and get back the checksum
  my $seq_obj = $exon->seq();
  my $seq_checksum;
  if($seq_obj) {
    $seq_checksum = $self->_insert_sequence($seq_obj->seq(), $session_pkg->{session_id});
  }

  my @loc_pieces = (
    $session_pkg->{assembly_id},$exon->seq_region_name(),
    $exon->seq_region_start(), $exon->seq_region_end(),
    $exon->seq_region_strand(),
  );

  my @loc_pieces_checksum = (
    $session_pkg->{assembly_name},$exon->seq_region_name(),
    $exon->seq_region_start(), $exon->seq_region_end(),
    $exon->seq_region_strand(),
  );

  my $utils = Bio::EnsEMBL::Tark::Utils->new();
  my $loc_checksum = $utils->checksum_array( @loc_pieces_checksum );
  my $exon_checksum = $utils->checksum_array( $loc_checksum, $seq_checksum );

  my $sth = $self->get_insert('exon');
  $sth->execute(
    $exon->stable_id(), $exon->version(), @loc_pieces, $loc_checksum,
    $exon_checksum, $seq_checksum, $session_pkg->{session_id}
  ) or $self->log->logdie("Error inserting exon: $DBI::errstr");
  my $exon_id = $sth->{mysql_insertid};

  # Apply tags to feature we've just inserted
  $tag->tag_feature($exon_id, 'exon');

  return ($exon_id, $exon_checksum);
} ## end sub _load_exon


=head2 _load_translation
  Arg [1]    : $translation : Bio::EnsEMBL::Translation
  Arg [2]    : $session_pkg : Bio::EnsEMBL::Tark::DB
  Arg [3]    : $tag         : Bio::EnsEMBL::Tark::Tag
  Description: For a translation generate and load the matching checksum
  Returntype :
  Exceptions : none
  Caller     : general

=cut

sub _load_translation {
  my ( $self, $translation, $session_pkg, $tag ) = @_;

  # Insert the sequence and get back the checksum
  my $seq_checksum = $self->_insert_sequence($translation->seq(), $session_pkg->{session_id});

  my @loc_pieces = (
    $session_pkg->{assembly_id}, $session_pkg->{transcript}->seq_region_name(),
    $translation->genomic_start(), $translation->genomic_end(),
    $session_pkg->{transcript}->seq_region_strand(),
  );

  my @loc_pieces_checksum = (
    $session_pkg->{assembly_name}, $session_pkg->{transcript}->seq_region_name(),
    $translation->genomic_start(), $translation->genomic_end(),
    $session_pkg->{transcript}->seq_region_strand(),
  );

  my $utils = Bio::EnsEMBL::Tark::Utils->new();
  my $loc_checksum = $utils->checksum_array( @loc_pieces_checksum );
  my $translation_checksum = $utils->checksum_array(
    $loc_checksum, $translation->stable_id(), $translation->version(), $seq_checksum
  );

  my $sth = $self->get_insert('translation');
  $sth->execute(
    $translation->stable_id(), $translation->version(), @loc_pieces, $loc_checksum,
    $translation_checksum, $seq_checksum, $session_pkg->{session_id}
  ) or $self->log->logdie("Error inserting translation: $DBI::errstr");
  my $translation_id = $sth->{mysql_insertid};

  $sth = $self->get_insert('translation_transcript');
  $sth->execute(
    $session_pkg->{transcript_id}, $translation_id, $session_pkg->{session_id}
  ) or $self->log->logdie("Error inserting translation_transcript: $DBI::errstr");

  # Apply tags to feature we've just inserted
  $tag->tag_feature($translation_id, 'translation');

  return;
} ## end sub _load_translation


=head2 _insert_sequence
  Arg [1]    : string - sequence
  Arg [2]    : integer - session_id
  Description: Load a sequence into the sequence tark schema table with the
               matching session_id and calculated checksum value
  Returntype : binary - sha1
  Exceptions : none
  Caller     : general

=cut

sub _insert_sequence {
  my ( $self, $sequence, $session_id ) = @_;

  my $utils = Bio::EnsEMBL::Tark::Utils->new();
  my $sha1 =  $utils->checksum_array( $sequence );

  my $sth = $self->get_insert('sequence');
  $sth->execute(
    $sha1, $sequence, $session_id
  ) or $self->log->logdie("Error inserting sequence: $DBI::errstr");

  return $sha1;
} ## end sub _insert_sequence


=head2 _fetch_name_id
  Description: Retrieve and formate the consortium assigned name_id.
  Returntype :
  Exceptions : none
  Caller     : general
=cut

sub _fetch_name_id {
  my ( $self, $gene, $consortium_name, $add_prefix ) = @_;

  foreach my $oxref (@{ $gene->get_all_object_xrefs() }) {
    if ( $oxref->dbname ne $consortium_name ) {
      next;
    }

    if ( defined $add_prefix and $add_prefix == 1 ) {
      return $consortium_name . q{:} . $oxref->primary_id;
    }

    return $oxref->primary_id;
  }

  return;
} ## end sub _fetch_name_id


=head2 genes_to_metadata_iterator
  Description: This is the place where you should try to get the gene iterator properly
  Returntype :
  Exceptions : none
  Caller     : general

=cut

sub genes_to_metadata_iterator {
  my ( $self, $dba, $source_name, $gene_ids ) = @_;
  my $ga = $dba->get_GeneAdaptor();

  if ( !defined $gene_ids ) {
    $gene_ids = $ga->_list_dbIDs('gene');
  }

  my $len          = scalar @{ $gene_ids };
  my $current_gene = 0;
  my $genes_i      = Bio::EnsEMBL::Utils::Iterator->new(
    sub {
      if ( $current_gene >= $len ) {
        return;
      }
      else {
        my $gene = $ga->fetch_by_dbID( $gene_ids->[ $current_gene++ ] );
        return $gene;
      }
    }
  );
  return $genes_i;
} ## end sub genes_to_metadata_iterator

1;
