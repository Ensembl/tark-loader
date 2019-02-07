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

package Bio::EnsEMBL::Tark::Schema::Result::Genome;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::EnsEMBL::Tark::Schema::Result::Genome

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

=head1 TABLE: C<genome>

=cut

__PACKAGE__->table("genome");

=head1 ACCESSORS

=head2 genome_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 1
  size: 128

=head2 tax_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 session_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_foreign_key: 1
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "genome_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "name",
  { data_type => "varchar", is_nullable => 1, size => 128 },
  "tax_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
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

=item * L</genome_id>

=back

=cut

__PACKAGE__->set_primary_key("genome_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<genome_idx>

=over 4

=item * L</name>

=item * L</tax_id>

=back

=cut

__PACKAGE__->add_unique_constraint("genome_idx", ["name", "tax_id"]);

=head1 RELATIONS

=head2 assemblies

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::Assembly>

=cut

__PACKAGE__->has_many(
  "assemblies",
  "Bio::EnsEMBL::Tark::Schema::Result::Assembly",
  { "foreign.genome_id" => "self.genome_id" },
  { cascade_copy => 0, cascade_delete => 0 },
);

=head2 assembly_aliases

Type: has_many

Related object: L<Bio::EnsEMBL::Tark::Schema::Result::AssemblyAlias>

=cut

__PACKAGE__->has_many(
  "assembly_aliases",
  "Bio::EnsEMBL::Tark::Schema::Result::AssemblyAlias",
  { "foreign.genome_id" => "self.genome_id" },
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
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:7U13wQ/V+hhT6UdVf4XKNQ

sub sqlt_deploy_hook {
  my ($self, $sqlt_table) = @_;

  $sqlt_table->add_index(name => 'fk_genome_1_idx', fields => ['session_id']);

  return;
}

# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
