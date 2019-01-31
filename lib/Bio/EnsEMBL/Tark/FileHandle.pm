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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Tark::FileHandle

=head1 SYNOPSIS

  use Bio::EnsEMBL::Tark::FileHandle;

  my $ft = Bio::FormatTranscriber->new(

=head1 DESCRIPTION



=cut

package Bio::EnsEMBL::Tark::FileHandle;

use strict;
use warnings;

use Carp;

use PerlIO::gzip;


=head2 get_file_handle
  Arg [1]    : $filename : string
  Description: Create a file handle for a given string. Capable of handling local
               files as well as URLs as well as gzipped files.
  Returntype : File handle
  Exceptions : none
  Caller     : general

=cut

sub get_file_handle {
  my $class = shift;
  my $filename = shift;

  my $mode = '<';

  if($filename =~ /^http/i) {
    $mode = '-|';
    $filename = 'curl -vs 2>/dev/null ' . $filename;
  }

  if($filename =~ /gz/i) {
    $mode .= ':gzip';
  }

  open my $fh, $mode, $filename or
    confess "Error opening file $filename: $!";

  return $fh;
} ## end sub get_file_handle

1;
