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

package Bio::EnsEMBL::Tark::Schema::Result::TranscriptReleaseTag;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::EnsEMBL::Tark::Schema::Result::TranscriptReleaseTag

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 COMPONENTS LOADED

=over 4

=item * L<DBIx::Class::InflateColumn::DateTime>

=back

=cut

__PACKAGE__->load_components("InflateColumn::DateTime");

=head1 TABLE: C<transcript_release_tag>

=cut

__PACKAGE__->table("transcript_release_tag");

=head1 ACCESSORS

=head2 transcript_release_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 feature_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 0

=head2 release_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 0

=head2 session_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "transcript_release_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "feature_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 0,
  },
  "release_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_foreign_key => 1,
    is_nullable => 0,
  },
  "session_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</transcript_release_id>

=back

=cut

__PACKAGE__->set_primary_key("transcript_release_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<transcript_release_idx>

=over 4

=item * L</feature_id>

=item * L</release_id>

=back

=cut

__PACKAGE__->add_unique_constraint("transcript_release_idx", ["feature_id", "release_id"]);

=head1 RELATIONS

=head2 feature

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Transcript>

=cut

__PACKAGE__->belongs_to(
  "feature",
  "Bio::EnsEMBL::Tark::Schema::Result::Transcript",
  { transcript_id => "feature_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);

=head2 release

Type: belongs_to

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::ReleaseSet>

=cut

__PACKAGE__->belongs_to(
  "release",
  "Bio::EnsEMBL::Tark::Schema::Result::ReleaseSet",
  { release_id => "release_id" },
  { is_deferrable => 1, on_delete => "NO ACTION", on_update => "NO ACTION" },
);


# Created by DBIx::Class::Schema::Loader v0.07045 @ 2019-01-22 14:34:33
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:hVBBsWX2gaXFGD8I9XPl0A

=head2 transcript_release_object

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::transcript_release_object>

=cut

__PACKAGE__->has_many(
  "transcript_release_object",
  "Bio::EnsEMBL::Tark::Schema::Result::TranscriptReleaseTagRelationship",
  { "foreign.transcript_release_object_id" => "self.transcript_release_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 transcript_release_subject

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::transcript_release_object>

=cut

__PACKAGE__->has_many(
  "transcript_release_subject",
  "Bio::EnsEMBL::Tark::Schema::Result::TranscriptReleaseTagRelationship",
  { "foreign.transcript_release_subject_id" => "self.transcript_release_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
