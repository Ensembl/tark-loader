use utf8;
package Bio::EnsEMBL::Tark::Schema::Result::OperonTranscript;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::EnsEMBL::Tark::Schema::Result::OperonTranscript

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

=head1 TABLE: C<operon_transcript>

=cut

__PACKAGE__->table("operon_transcript");

=head1 ACCESSORS

=head2 operon_transcript_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 1
  size: 64

=head2 stable_id_version

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 operon_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=head2 transcript_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=head2 session_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "operon_transcript_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "stable_id",
  { data_type => "varchar", is_nullable => 1, size => 64 },
  "stable_id_version",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "operon_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
  "transcript_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 1,
  },
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

=item * L</operon_transcript_id>

=back

=cut

__PACKAGE__->set_primary_key("operon_transcript_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<operon_transcript_idx>

=over 4

=item * L</operon_id>

=item * L</transcript_id>

=back

=cut

__PACKAGE__->add_unique_constraint("operon_transcript_idx", ["operon_id", "transcript_id"]);

=head1 RELATIONS

=head2 operon

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Operon>

=cut

__PACKAGE__->belongs_to(
  "operon",
  "Bio::EnsEMBL::Tark::Schema::Result::Operon",
  { operon_id => "operon_id" },
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

=head2 transcript

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Transcript>

=cut

__PACKAGE__->belongs_to(
  "transcript",
  "Bio::EnsEMBL::Tark::Schema::Result::Transcript",
  { transcript_id => "transcript_id" },
  {
    is_deferrable => 1,
    join_type     => "LEFT",
    on_delete     => "NO ACTION",
    on_update     => "NO ACTION",
  },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2019-01-22 14:34:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:OQMS9DmOJmFV51h4Ef8BRA

sub sqlt_deploy_hook {
  my ($self, $sqlt_table) = @_;

  $sqlt_table->add_index(name => 'stable_id', fields => ['stable_id', 'stable_id_version']);
}

# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
