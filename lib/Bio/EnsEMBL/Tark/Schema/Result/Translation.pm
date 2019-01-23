use utf8;
package Bio::EnsEMBL::Tark::Schema::Result::Translation;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::EnsEMBL::Tark::Schema::Result::Translation

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<DBIx::Class::InflateColumn::DateTime>

=back

=cut

__PACKAGE__->load_components("InflateColumn::DateTime");

=head1 TABLE: C<translation>

=cut

__PACKAGE__->table("translation");

=head1 ACCESSORS

=head2 translation_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 0
  size: 64

=head2 stable_id_version

  data_type: 'tinyint'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 assembly_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=head2 loc_start

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 loc_end

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 loc_strand

  data_type: 'tinyint'
  is_nullable: 1

=head2 loc_region

  data_type: 'varchar'
  is_nullable: 1
  size: 42

=head2 loc_checksum

  data_type: 'binary'
  is_nullable: 1
  size: 20

=head2 translation_checksum

  data_type: 'binary'
  is_nullable: 1
  size: 20

=head2 seq_checksum

  data_type: 'binary'
  is_foreign_key: 1
  is_nullable: 1
  size: 20

=head2 session_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "translation_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "stable_id",
  { data_type => "varchar", is_nullable => 0, size => 64 },
  "stable_id_version",
  { data_type => "tinyint", extra => { unsigned => 1 }, is_nullable => 0 },
  "assembly_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
  "loc_start",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "loc_end",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "loc_strand",
  { data_type => "tinyint", is_nullable => 1 },
  "loc_region",
  { data_type => "varchar", is_nullable => 1, size => 42 },
  "loc_checksum",
  { data_type => "binary", is_nullable => 1, size => 20 },
  "translation_checksum",
  { data_type => "binary", is_nullable => 1, size => 20 },
  "seq_checksum",
  { data_type => "binary", is_foreign_key => 1, is_nullable => 1, size => 20 },
  "session_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</translation_id>

=back

=cut

__PACKAGE__->set_primary_key("translation_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<translation_chk>

=over 4

=item * L</translation_checksum>

=back

=cut

__PACKAGE__->add_unique_constraint("translation_chk", ["translation_checksum"]);

=head1 RELATIONS

=head2 assembly

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Assembly>

=cut

__PACKAGE__->belongs_to(
  "assembly",
  "Bio::EnsEMBL::Tark::Schema::Result::Assembly",
  { assembly_id => "assembly_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);

=head2 seq_checksum

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Sequence>

=cut

__PACKAGE__->belongs_to(
  "seq_checksum",
  "Bio::EnsEMBL::Tark::Schema::Result::Sequence",
  { seq_checksum => "seq_checksum" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);

=head2 session

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Session>

=cut

__PACKAGE__->belongs_to(
  "session",
  "Bio::EnsEMBL::Tark::Schema::Result::Session",
  { session_id => "session_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);

=head2 translation_release_tags

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::TranslationReleaseTag>

=cut

__PACKAGE__->has_many(
  "translation_release_tags",
  "Bio::EnsEMBL::Tark::Schema::Result::TranslationReleaseTag",
  { "foreign.feature_id" => "self.translation_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 translation_transcripts

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::TranslationTranscript>

=cut

__PACKAGE__->has_many(
  "translation_transcripts",
  "Bio::EnsEMBL::Tark::Schema::Result::TranslationTranscript",
  { "foreign.translation_id" => "self.translation_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2019-01-22 14:34:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Ez3MGYP6ig+907cGcF8GdA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
