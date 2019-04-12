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

package Bio::EnsEMBL::Tark::Hive::PipeConfig::TarkSQL;

use Moose;


=head1 _gene_grouping_template_SQL
  Description: Default SQL for the gene grouping queries. This is an internal
               function and should only be accessed through the gene_grouping,
               gene_grouping_exclusion and gene_grouping_inclusion that will
               handle the substitution of the WHERE clauses.
=cut

sub _feature_template_SQL {
  my ($self) = @_;

  my $sql = (<<'SQL');
    SELECT
      COUNT(*)
    FROM
      #FEATURE#
SQL

  return $sql;
} ## end sub _gene_grouping_template_SQL


=head2 feature_count
  Description: Query for getting lists of all genes split into n batches. Genes
               are randomly assigned to each group.
=cut

sub feature_count {
  my ( $self, $feature ) = @_;

  my $sql = $self->_feature_template_SQL();

  $sql =~ s/#FEATURE#/$feature/g;

  return $sql;
} ## end sub feature_count


=head1 _feature_release_template_SQL
  Description: Default SQL for the gene grouping queries. This is an internal
               function and should only be accessed through the gene_grouping,
               gene_grouping_exclusion and gene_grouping_inclusion that will
               handle the substitution of the WHERE clauses.
=cut

sub _feature_release_template_SQL {
  my ($self) = @_;

  my $sql = (<<'SQL');
    SELECT
      #FEATURE#_release_tag.release_id,
      COUNT(*)
    FROM
      #FEATURE#_release_tag
      JOIN release_set ON #FEATURE#_release_tag.release_id=release_set.release_id
      JOIN release_source ON #FEATURE#_release_tag.source_id=release_source.source_id
    WHERE
      release_set.shortname = ? AND
      release_source.short_name = ?
    GROUP BY
      #FEATURE#_release_tag.release_id
SQL

  return $sql;
} ## end sub _feature_release_template_SQL


=head2 feature_release_count
  Arg [1]    : $feature - string
  Description: Query for getting counts of all features
=cut

sub feature_release_count {
  my ( $self, $feature ) = @_;

  my $sql = $self->_feature_release_template_SQL();

  $sql =~ s/#FEATURE#/$feature/g;

  return $sql;
} ## end sub feature_release_count


=head1 _feature_release_template_SQL
  Description: Default SQL for the gene grouping queries. This is an internal
               function and should only be accessed through the gene_grouping,
               gene_grouping_exclusion and gene_grouping_inclusion that will
               handle the substitution of the WHERE clauses.
=cut

sub _feature_diff_template_SQL {
  my ($self) = @_;

  my $sql = (<<'SQL');
    SELECT
      COUNT(*)
    FROM
      (
        SELECT
          #FEATURE#.stable_id,
          #FEATURE#.stable_id_version,
          f_tag.feature_id,
          rst.shortname,
          rst.description,
          rst.assembly_id
        FROM
          #FEATURE#
          JOIN #FEATURE#_release_tag AS f_tag ON (#FEATURE#.#FEATURE#_id=f_tag.feature_id)
          JOIN release_set AS rst ON (f_tag.release_id=rst.release_id)
          JOIN release_source AS rso ON (rst.source_id=rso.source_id)
        WHERE
          rst.shortname=? AND
          rso.shortname=?
      ) AS v0
      #DIRECTION# JOIN (
        SELECT
          #FEATURE#.stable_id,
          #FEATURE#.stable_id_version,
          f_tag.feature_id,
          rst.shortname,
          rst.description,
          rst.assembly_id
        FROM
          #FEATURE#
          JOIN #FEATURE#_release_tag AS f_tag ON (#FEATURE#.#FEATURE#_id=f_tag.feature_id)
          JOIN release_set AS rst ON (f_tag.release_id=rst.release_id)
          JOIN release_source AS rso ON (rst.source_id=rso.source_id)
        WHERE
          rst.shortname=? AND
          rso.shortname=?
      ) AS v1 ON (v0.stable_id=v1.stable_id)
    WHERE
      #SET#
SQL

  return $sql;
} ## end sub _feature_release_template_SQL


=head2 feature_diff_count
  Arg [1]    : $feature - string
  Arg [2]    : $direction - string (removed|gained)
  Description: Query for getting counts of all features
=cut

sub feature_diff_count {
  my ( $self, $feature, $direction ) = @_;

  my $sql = $self->_feature_diff_template_SQL();

  $sql =~ s/#FEATURE#/$feature/g;

  if ( $direction eq 'removed' ) {
    $sql =~ s/#DIRECTION#/LEFT/g;
    $sql =~ s/#SET#/v1.stable_id IS NULL/g;
  } elsif ($direction eq 'gained') {
    $sql =~ s/#DIRECTION#/RIGHT/g;
    $sql =~ s/#SET#/v0.stable_id IS NULL/g;
  } else {
    $sql =~ s/#DIRECTION#//g;
    $sql =~ s/#SET#/v0.stable_id_version!=v1.stable_id_version/g;
  }

  return $sql;
} ## end sub feature_diff_count


=head2 insert_stats
  Arg [1]    : $feature - string
  Arg [2]    : $direction - string (removed|gained)
  Description: Query for getting counts of all features
=cut

sub insert_stats {
  my ( $self, $json ) = @_;

  my $sql = (<<'SQL');
    INSERT INTO release_stats (release_id, json)
    SELECT
      release_set.release_id,
      '#JSON#' AS json
    FROM
      release_set
      JOIN release_source ON release_set.source_id=release_source.source_id
    WHERE
      release_source.shortname=? AND
      release_set.shortname=?
SQL

  $sql =~ s/#JSON#/$json/g;

  return $sql;
} ## end sub feature_diff_count

1;
