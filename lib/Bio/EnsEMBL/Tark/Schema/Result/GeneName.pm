use utf8;
package Bio::EnsEMBL::Tark::Schema::Result::GeneName;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::EnsEMBL::Tark::Schema::Result::GeneName

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

=head1 TABLE: C<gene_names>

=cut

__PACKAGE__->table("gene_names");

=head1 ACCESSORS

=head2 gene_names_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 external_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 name

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 source

  data_type: 'varchar'
  is_nullable: 1
  size: 32

=head2 primary_id

  data_type: 'tinyint'
  is_nullable: 1

=head2 session_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "gene_names_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "external_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "name",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "source",
  { data_type => "varchar", is_nullable => 1, size => 32 },
  "primary_id",
  { data_type => "tinyint", is_nullable => 1 },
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

=item * L</gene_names_id>

=back

=cut

__PACKAGE__->set_primary_key("gene_names_id");

=head1 RELATIONS

=head2 genes

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Gene>

=cut

__PACKAGE__->has_many(
  "genes",
  "Bio::EnsEMBL::Tark::Schema::Result::Gene",
  { "foreign.hgnc_id" => "self.external_id" },
  { cascade_copy => 0, cascade_delete => 0 },
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


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2019-01-22 14:34:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:lXDy1JKemCgTXGtWz5xRbw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
