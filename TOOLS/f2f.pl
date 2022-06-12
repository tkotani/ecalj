#!/usr/bin/perl -w
# Author: Colby Lemon (colbylemon@gmail.com)
# "Helpful" refactoring provided by Bob Apthorpe (apthorpe+cpan@cynistar.net)
# Version 0.95

# This program converts Fortran 77 code to Fortran 90
# USAGE: f2f [inputfile [outputfile]]

# If no arguments are given, standard input and output are assumed
# It will do the following:
# 1) Convert comments to f90 style comments
# 2) Change line continuation character to '&'
# 3) Change relational operators to ==, >=, etc.
# 4) Convert variable declarations to f90 :: style (see *** below)
# 5) Terminate do loops with 'end do'
# 6) Convert some 'go to' statements to 'cycle'
# 7) Append subroutine names to 'end' statements
# 8) Indent the resulting code

# *** WATCH OUT FOR THIS ***:
# Double Precision real variables are considered redundant in
# fortran 90, so this code tries to replace them.  Double precision
# variables are "kind=8" on most compilers, but this is not a standard.
# You may need to edit the $k assignment below so that $k=2 if you use,
# e.g., the NAG compiler

# For simplicity, this code sometimes assumes we are dealing with
# a "pure" Fortran 77 input file.  However, many people use mixed
# Fortran 77/90/95 code.  This may or may not cause problems.
# Give it a try and see

# TODO Add command line options for certain preferences (KIND value)
# TODO Handle non-standard F77 code better
# TODO Deal with END DO statements preceded by a numeric label
# TODO Remind people to convert their INCLUDEd code? [Done. Sort of.]
# TODO Ensure multiple ENDDO statements for nested DO statements?
# TODO Perform better on mixed Fortran77/90 code
# TODO Automatically convert to UNIX line endings, or allow Windows line endings (can perl somehow do this automaticall?)
# TODO Collect a library of fortran 77 code, to do a regular compile and test on, on multiple computers.
# TODO Think about indentation and any other user preferences, and give command line options
# TODO Update the copyright?

#    Copyright 2003 Colby Lemon
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#    See LICENSE.txt which should have been shipped with this program.
#
use warnings;
use strict;

use Carp;
use Data::Dumper;
use English qw(-no_match_vars);

# use File::Slurp;
use Getopt::Long;
use IO::File;
use Pod::Usage;

our $CODENAME = q{f2f};
our $VERSION  = q{0.95};

# TODO Set these options as defaults to be overridden by command line options
# TODO Read options from command line
#
#     Note that 'dp-to-kind' can be text (a parameter name) for
# parameterizing precision (see MKKIND.F90 at
# http://www.sesp.cse.clrc.ac.uk/Publications/Legacy-Software/Legacy-Software/node38.html)
#
#     $rh_opt contains options such as tab stops, base indentation level, &c.
#
#     Note that although DOUBLE PRECISION is deprecated, REAL*8 is not
# even standard. REAL*8 is a Digital (VAX) Fortran extension.
#
my $rh_default = {

    # Convert DOUBLE PRECISION to REAL*<dp-to-star-kind>
    'dp-to-star-kind' => undef,

    # Convert DOUBLE PRECISION to REAL (KIND=<dp-to-kind>)
    'dp-to-kind' => undef,

    # Base indent is one tab stop
    'base-indent' => 1,

    # Tab stop is four spaces
    'tab' => 4,
};

# Set options to defaults - this is in preparation for handling a
# config file or command-line options
my $rh_opt = {};
%{$rh_opt} = %{$rh_default};

process_options($rh_opt,);

my $status = main($rh_opt);

exit $status;

sub main {
    my $rh_opt = shift;
    my $status = 0;

    # reference to list of source code lines
    my $rl_src = [];
    
    if ($rh_opt->{verbose}) {
        print qq{Processing $ARGV[0]\n};
    }

    # slurp source code in from file or STDIN
    $status += read_source($rl_src, $ARGV[0]);

    # Free the source from fixed form
    $status += liberate($rh_opt, $rl_src);

    # Now make changes to the code
    $status += convert($rh_opt, $rl_src);

    # Indent the code
    $status += indent($rh_opt, $rl_src);

    $status += write_source($rh_opt, $rl_src);

    return $status;
}

# =-=-=-=-=-=-=-=-=-=-=-= End of Main Program =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

sub process_options {
    my $rh_opt = shift;

    ### set option types and arguments

    my @optlst = (
        'tab=i',
        'base-indent=i',
        'dp-to-star-kind=i',
        'dp-to-kind=s',
        'suffix-includes=s',
        'help',
        'usage',
        'version',
        'verbose',
    );

    ### map options to variables

    for my $kk (@optlst) {
        my ($nkk, $junk) = split(m/[=!]/, $kk, 2);
    }

    ## Parse options
    GetOptions($rh_opt, @optlst) || pod2usage(-verbose => 0);

    ### print usage, version, or ownership messages.

    if (   $rh_opt->{'help'}
        || $rh_opt->{'usage'}
        || $rh_opt->{'version'})
    {

        if ($rh_opt->{'help'}) {
            pod2usage(-verbose => 1);
        }

        if ($rh_opt->{'usage'}) {
            pod2usage(-verbose => 0);
        }

        if ($rh_opt->{'version'}) {
            print STDERR qq{$CODENAME, version $VERSION\n};
        }

        exit 0;
    }

    return;
}

sub read_source {
    my $rl_src = shift;
    my $ifn    = shift;

    my $status = 0;

    my $ifh;
    if ($ifn) {

        # if a filename is specified
        if (-f $ifn) {

            # open it for reading
            $ifh = IO::File->new($ifn, q{r});
            if ($ifh) {

                # read it all into a buffer
                @{$rl_src} = <$ifh>;
                $ifh->close();
            } else {
                $status = -1;
                carp(qq{ERROR: Cannot read file '$ifn' - $OS_ERROR\n});
                exit $status;
            }
        } else {
            $status = -2;
            croak(qq{ERROR: File '$ifn' does not exist.\n});
            exit $status;
        }
    } else {

        # otherwise read from standard input
        @{$rl_src} = <>;
    }
    return $status;
}

sub write_source {
    my $rh_opt = shift;
    my $rl_src = shift;

    my $status = 0;

    # print the Fortran 90 code
    my $ofh;

    if ($ARGV[1]) {

# if a second argument is specified, it is the filename of the newly-reformatted code
        my $ofn = $ARGV[1];

        # open it for output
        $ofh = IO::File->new($ofn, q{w});

        if ($ofh) {

            # Send to file
            print $ofh @{$rl_src};
            $ofh->close();
        } else {
            print STDERR
                qq{Warning: Cannot write file to '$ofn' - $OS_ERROR\n};
            $status = -4;

            # Consolation prize - send it to STDOUT
            print @{$rl_src};
        }
    } else {

        # send it to STDOUT
        print @{$rl_src};
    }

    return $status;
}

sub liberate {
    my $rh_opt   = shift;
    my $rl_src   = shift;
    my $last_idx = scalar(@{$rl_src}) - 1;

    my $status = 0;

libloop:
    for my $lineno (0 .. $last_idx) {
        my $card = $rl_src->[$lineno];

        # throw out trailing whitespace
        $card =~ s/[ \t]+$//;

        # skip blank lines
        next libloop if ($card =~ m/^\n/);

        # Convert comments
        # change 'c...' to '!...'
        $card =~ s/^[c*]/!/io;

        # use '&' for line continuation

        # if character in column 6...
        if ($card =~ s/^ {5}[^0 \t]//o) {

            # find last statement

            # (skip blank lines and comments)
            my $j = 1;
            while ($rl_src->[ $lineno - $j ] =~ /^[\n!]/) {
                ++$j;
            }

            $rl_src->[ $lineno - $j ] =
                add_continuation_marker($rl_src->[ $lineno - $j ]);
        }

        # Get rid of fixed source formatting spaces

        # throw out leading spaces
        $card =~ s/^\s+//;

        # compress spaces after numeric labels
        $card =~ s/^([0-9]+)\s{2,}(.+)/$1 $2/;

        # put modified line back into source list
        $rl_src->[$lineno] = $card;
    }

    return $status;
}

sub convert {
    my $rh_opt   = shift;
    my $rl_src   = shift;
    my $last_idx = scalar(@{$rl_src}) - 1;

    my $status = 0;

    my @newsrc = ();

    # TODO Refactor away implied variable
    @_ = @{$rl_src};

    my $subname;
    my $cn;
    my $linect = 0;

    # initialize array to store loop labels
    my @label = (0);

    my $sfx;
    if ($rh_opt->{'suffix-includes'}) {
        $sfx = $rh_opt->{'suffix-includes'};
    }    

cloop:
    for (my $lineno = 0; $_ = $_[$lineno]; $_[ $lineno++ ] = $_) {

        #        print STDERR q{convert: Line } . $linect . qq{ =? $lineno\n};
        $linect++;

        # skip blank lines and comments
        next cloop if /^(\n|!)/;

    cblock:
        {
            # add suffix to includes files 
            if ($sfx && /^\s*include\b/i) {
                s/([a-z][a-z0-9_.]*)(['"]\s*$)/${1}${sfx}${2}/i;
                print STDERR qq{Note: Remember to convert $1 to $1$sfx\n};
                last cblock;
            }

            # replace .eq. , .gt. , etc. with ==, >, etc., and add spaces
            if (/^[0-9]*\s*if/i || /^else\s?if/i || $cn) {
                s/ ?\.\s*lt\s*\. ?/ < /ig;
                s/ ?\.\s*eq\s*\. ?/ == /ig;
                s/ ?\.\s*gt\s*\. ?/ > /ig;
                s/ ?\.\s*le\s*\. ?/ <= /ig;
                s/ ?\.\s*ne\s*\. ?/ \/= /ig;
                s/ ?\.\s*ge\s*\. ?/ >= /ig;

                # add spaces and capitalize logical operators
                s/ ?\.\s*and\s*\. ?/ .AND. /ig;
                s/ ?\.\s*or\s*\. ?/ .OR. /ig;
                s/ ?\.\s*not\s*\. ?/ .NOT. /ig;
                s/ ?\.\s*true\s*\. ?/ .TRUE. /ig;
                s/ ?\.\s*false\s*\. ?/ .FALSE. /ig;
                $cn = /&$/;

                last cblock;
            }

            # add program unit names to end statements
            # this includes subroutine, function, module, and block data
            if ((/^(
                    (subroutine|module|block data)
                    [ ][\w_\d]*)/ix)  # subroutine, module, block data
                || (/^\w*\s*(function [\w_\d]*)/i)    # function
                )
            {

                # save program unit name
                $subname = $1;
                last cblock;
            }

            # double precision function -> real*8 function
            if ($rh_opt->{'dp-to-star-kind'}) {
                my $dpkind = $rh_opt->{'dp-to-star-kind'};

                # TODO Adjust case
                if (s/^double precision (function [\w_\d]*)/real*$dpkind $1/i)
                {
                    $subname = $1;
                    last cblock;
                }
            } elsif ($rh_opt->{'dp-to-kind'}) {
                my $dpkindspec = $rh_opt->{'dp-to-kind'};

                # TODO Adjust case
                if (s/^double precision (function [\w_\d]*)/real (kind=$dpkindspec) $1/i
                    )
                {
                    $subname = $1;
                    last cblock;
                }
            }

            if (/^end$/i) {    # if end of program unit
                if ($subname) {
                    my $endcmd = 'END';
                    if ($subname =~ m/[a-z]/o) {
                        $endcmd = 'end';
                    }
                    s/^end$/$endcmd $subname/i;    # append subroutine name

                    $subname = '';
                } else {

                    # TODO Adjust case
                    s/^end$/END PROGRAM/i;
                }
                last cblock;
            }

            # Update declarations

            # Add :: to character variables
            if (s/^(character[ \t]*\*[1-9][0-9]*) /$1 :: /i) {

                # *123 -> (123)
                s/\*([1-9][0-9]*)/($1)/g;

                # TODO Refactor away enmanglement of the loop invariant. Ow.
                # Keep going if line continues
                while (/&$/) {
                    $_[$lineno] = $_;
                    $_ = $_[ ++$lineno ];

                    # *123 -> (123)
                    s/\*([0-9]*)/($1)/g;
                }
                last cblock;
            }

            {

                # add :: to real and integer
                # (as long as there's no : immediately downstream...)
                s/^(real|integer) (?![:]|FUNCTION)/$1 :: /i
                    && next cloop;

                # add :: to real*4
                s/^(real\*[1-9][0-9]*) /$1 :: /i
                    && next cloop;

                # add :: to integer*4
                s/^(integer\*[1-9][0-9]*) /$1 :: /i
                    && next cloop;

                # add :: to character*4
                s/^(character\*([1-9][0-9]*)) /$1 (LEN=$2) :: /i
                    && next cloop;

                # add :: to double precision
                s/^(double precision)/$1 ::/i;

                if ($rh_opt->{'dp-to-star-kind'}) {

                    # Replace DOUBLE PRECISION with REAL*8
                    # TODO Adjust case
                    my $dpkind = $rh_opt->{'dp-to-star-kind'};
                    s/^double precision(\s*::)?(?!.*[:])/REAL*$dpkind ::/i;
                } elsif ($rh_opt->{'dp-to-kind'}) {

                    # Replace DOUBLE PRECISION with REAL (KIND=_kind_)
                    # TODO Adjust case
                    my $dpkindspec = $rh_opt->{'dp-to-kind'};
                    s/^double precision(\s*::)?(?!.*[:])/REAL (KIND=$dpkindspec) ::/i;
                }

                # add :: to logical and complex
                s/^(logical|complex) (?!.*[:]|FUNCTION)/$1 :: /i;

                last cblock;
            }
        }    # cblock

        # change do loop terminator to 'end do'
        if (/^[0-9]*[ \t]*do[ \t]+([0-9]+)/i) {    # do loop
            push @label, $1;                       # save the label

        } elsif (/^([0-9]+)[ \t]/i && ($1 eq $label[-1])) {    # end of loop
            my $label = $1;    # save label number
            unless (s/continue$/END DO/i) {    # if not a continue statement,
                s/^([0-9]+)[ \t]+//;           # remove the label
                while (/&$/) {                 # while line continues
                    $_[ $lineno++ ] = $_;
                    $_ = $_[$lineno];    # go to the next line
                }

             # add 'END DO' line (clearer than ending loop at final statement)
             # TODO Replace with array splice()
                @_ = (
                    @_[ 0 .. $lineno ],
                    qq{$label END DO\n},
                    @_[ $lineno + 1 .. $#_ ]);
            }
            do {
                pop @label;    # remove labels from @label
                } while ($label[-1] == $label)

                #              } while ( $label[$#label] == $label )
                ;              # that belong to this do loop

            # replace some 'go to' statements with 'cycle'
            push @_, $lineno;    # store line number
            do {
                --$lineno;
            } while ($_[$lineno] =~ /^!/ || $_[$lineno] =~ /^\n/);

           # if line before end of loop is 'continue', search loop for 'go to'
            if ($_[$lineno] =~ s/^([0-9]+)\s+continue//i) {
                $label = $1;
                until ($_[ --$lineno ] =~ /^do[ \t]/i) {

                    # replace with 'cycle'
                    $_[$lineno] =~ s/go\s*to\s+$label$/cycle/i;

                }
            }

            # restore line number
            $lineno = pop @_;
        }
    }    # end of for loop

    # TODO Refactor away implied variable
    @{$rl_src} = @_;

    return $status;
}    # end of sub convert

sub indent {
    my $rh_opt = shift;
    my $rl_src = shift;

    my $status = 0;

    my $tab = $rh_opt->{tab};

    my $ccn;
    my $cnlabel;
    my @clabel;
    my $j = $rh_opt->{'base-indent'};

    my $linect = 0;

    my @newsrc = ();

iloop:
    foreach my $card (@{$rl_src}) {
        $linect++;
        if ($card =~ m/^\n/) {
            push @newsrc, $card;
            next iloop;
        }

        # make comments look nice
    iblock:
        {
            if ($card =~ m/^!/) {

                # insert space in worded comments
                # remove comments from blank lines
                $card =~ s/^![^\s]([\d\w]+)/! $1/ || $card =~ s/^!$//;

                # indent the line
                $card = q{ } x ($tab * ($j - 1)) . $card;
                last iblock;
            }

            # if loop indentation
            if ($card =~ m/^[0-9]* *if *\(/i || $ccn) {

                # single or continued-line if?
                $ccn = ($card =~ m/&$/);

                if ($card =~ m/\)\s*then/i) {

                    # indent line, increment $j
                    $card = q{ } x ($tab * $j) . $card;
                    $j++;
                } elsif ($card =~ /\bthen\b/i) {

                    # indent line, increment $j
                    $card = q{ } x ($tab * $j) . $card;
                    $j++;
                } else {

                    # indent line
                    $card = q{ } x ($tab * $j) . $card;
                }
                last iblock;
            }

            # END IF - close IF block
            if ($card =~ m/^[0-9]* *end *if/i) {

                # decrement $i, indent line
                $j--;
                $card = q{ } x ($tab * $j) . $card;
                last iblock;
            }

            # Outdent ELSE
            if ($card =~ /^[0-9]* *else/i) {

                # indent else with one less tab
                $card = q{ } x ($tab * ($j - 1)) . $card;
                last iblock;
            }

            # DO loop indentation (numeric label)
            if ($card =~ /^[0-9]*\s*do[ \t]+([0-9]*)/i) {
                if ($1) {

                    # save the label
                    push @clabel, $1;
                }

                # indent line, increment $i
                $card = q{ } x ($tab * $j) . $card;
                $j++;
                last iblock;
            }

            # DO loop indentation (text label)
            if ($card =~ /^[0-9]*\s*([a-z][^: ]*):\s*do\s/i) {
                $cnlabel = $1;
                if ($cnlabel) {
                    chomp $cnlabel;

                    # save the label
                    push @clabel, $1;

## Diagnostic
#                    print STDERR q{+}
#                        . ($j + 1)
#                        . qq{ - Text loop label $1 found\n};

                    # indent line, increment $i
                    $card = q{ } x ($tab * $j) . $card;
                    $j++;
                }
                last iblock;
            }

            # Labeled 'END DO' statement (numeric label)
            if ($card =~ m/^([0-9]+)[ \t]end *do$/i) {

                # save label
                $cnlabel = $1;

                # decrement $j as needed
            labelscan:
                while ($cnlabel eq $clabel[-1]) {
                    my $clabelct = $#clabel;
                    pop @clabel;
                    $j--;
                    last labelscan if ($#clabel < 0);
                }

                # continue -> end do
                $card = q{ } x ($tab * $j) . $card;
                last iblock;
            }

            # Labeled 'END DO' statement (text label)
            if ($card =~ m/^\s*end *do *([a-z]\w*)$/i) {

                # save label
                $cnlabel = $1;
                chomp $cnlabel;

                # decrement $j as needed
                if ($#clabel >= 0) {
                labelscan:
                    while ($cnlabel eq $clabel[-1]) {
                        my $clabelct = $#clabel;
                        pop @clabel;
                        $j--;
                        last labelscan if ($#clabel < 0);
                    }
                }

## Diagnostic
#                print STDERR q{-}
#                    . ($j + 1)
#                    . qq{ - Text loop label $cnlabel found\n};

                # indent
                $card = q{ } x ($tab * $j) . $card;
                last iblock;
            }

            # Bare 'END DO' statement
            if ($card =~ m/^end *do/i) {

                # decrement $j, indent line
                $j--;
                $card = q{ } x ($tab * $j) . $card;
                last iblock;
            }

            # WHERE construct indentation
            if ($card =~ m/^[0-9]*\s*where[ \t]+([0-9]*)/i)
            {    # where construct
                    # save the label
                if ($1) {
                    push @clabel, $1;
                }

                # indent line & increment $i
                $card = q{ } x ($tab * $j) . $card;
                $j++;
                last iblock;
            }

            # END WHERE
            if ($card =~ m/^end *where/i) {

                # decrement $j and indent line
                $j--;
                $card = q{ } x ($tab * $j) . $card;
                last iblock;
            }

            # default behavior - indent the line
            $card = q{ } x ($tab * $j) . $card;
        }

        # push modified card onto return stack
        push @newsrc, $card;

        if ($j < $rh_opt->{'base-indent'}) {
            print STDERR 'Warning: Indent level $j < $rh_opt->{base-indent} ('
                . $j . ' < '
                . $rh_opt->{'base-indent'} . ').'
                . '  There could be indentation problems at line '
                . $linect . qq{\n};
            $j      = $rh_opt->{'base-indent'};
            $status = -8;
        }
    }

    # replace original code with modified code
    @{$rl_src} = @newsrc;

    return $status;
}

# Finding the end of a line containing a string literal with a trailing
# comment is a somwehat complex problem, especially considering there
# are two different string delimiters. Happily, it's not as complex as
# finding matching parentheses. This subroutine implements a state
# machine to determine if a comment marker (!) is outside a string
# literal, thus ending a line and indicating where a trailing 
# continuation marker should be inserted.
sub add_continuation_marker {
    my $line = shift;

    chomp $line;
    
    # transitions to state 0 that may indicate a problem
    my %warnable = (
        2 => 1,
        3 => 1,
        5 => 1,
        6 => 1,
    );
    
    # state machine transition table
    # key 1: initial state
    # key 2: search pattern
    # value: new state
    my %t;
    $t{1}{qr{'}} = 2;
    $t{1}{qr{"}} = 5;
    $t{1}{qr{!}} = 0;
    $t{1}{qr{[^'"!]}} = 1;
    $t{2}{qr{'}} = 1;
    $t{2}{qr{[^']}} = 3;
    $t{3}{qr{'}} = 4;
    $t{3}{qr{[^']}} = 3;
    $t{4}{qr{'}} = 3;
    $t{4}{qr{!}} = 0;
    $t{4}{qr{[^'!]}} = 1;
    $t{5}{qr{"}} = 1;
    $t{5}{qr{[^"]}} = 6;
    $t{6}{qr{"}} = 7;
    $t{6}{qr{[^"]}} = 6;
    $t{7}{qr{"}} = 6;
    $t{7}{qr{!}} = 0;
    $t{7}{qr{[^'!]}} = 1;
    
    my $lastpos = (length $line) - 1;
    my $state = 1;
    my $i = 0;
    my $c = substr $line, $i, 1;

    my $eol = 0;
    stateloop:
    while ($state > 0) {
        if ($i > $lastpos) {
            $eol = 1;
            $state = 0;
            if (exists $warnable{$state}) {
                carp qq{Warning: Premature end of continued undelimited character string:\n$line\n}; 
            }
            last stateloop;
        }

        my $found = 0;
        scanloop:
        for my $tpat (keys %{$t{$state}}) {
            if ($c =~ m/${tpat}/) {
                $state = $t{$state}{$tpat};
                $i++;
                $c = substr $line, $i, 1;
                $found = 1;
                last scanloop;
            }
        }

        if ($found == 0) {
                # We should never get here
                carp qq{Warning: Transition not found (state = $state, i = $i, c = $c):\n$line\n};
                # bail
                $eol = 1;
                last stateloop;
        }
    }
    
    if ($eol) {
        $line .= q{ &};
    } else {
        $line = substr($line, 0, $i) . q{ & } . substr($line, $i) 
    }
    $line .= qq{\n};
    
    #TODO Warn if continued line is too long
    
    return $line;
}

__END__
# Add pod here
=pod

=head1 NAME

    f2f - FORTRAN-to-Fortran converter

=head1 SYNOPSIS

    f2f SOURCE.FOR dest.f90

=head1 VERSION

    $Revision: 0.95 $

=head1 DESCRIPTION

This program converts FORTRAN 77 code to Fortran 90. Mostly.

It will do the following:

  1) Convert comments to f90 style comments (! anywhere vs C in column 1)
  2) Change line continuation character to '&' at EOL of lines to continue (vs mark in column 6 in lines which continue)
  3) Change relational operators to from .EQ. to ==, .GE. to >=, etc.
  4) Convert variable declarations to f90 :: style (see BUGS AND LIMITATIONS below)
  5) Terminate do loops with 'end do'
  6) Convert some 'go to' statements to 'cycle'
  7) Append subroutine and function names to 'end' statements
  8) Indent the resulting code

If no arguments are given, standard input and output are assumed

=head1 USAGE

 Usage: f2f [source [destination]]
     --tab <n>               Converts tabs to <n> spaces, indents
                             <n> spaces per level (default: 4)
     --base-indent <n>       Sets initial indent level to <n>
                             (default: 1)
     --dp-to-star-kind <n>   Convert DOUBLE PRECISION to 
                             REAL*<dp-to-star-kind>
                             (no default)
     --dp-to-kind <str>      Convert DOUBLE PRECISION to 
                             REAL (KIND=<dp-to-kind>)
                             (no default)
     --suffix-includes <str> Adds suffix <str> to file names listed
                             in INCLUDE statements (no default)
     --verbose               Show information about code progress
     --help                  Displays this message
     --usage                 Displays program usage statement
     --version               Shows version info

=head1 REQUIRED ARGUMENTS

None.

=head1 OPTIONS

=over 4

=item --tab <n>

Converts tabs to <n> spaces, indents <n> spaces per level (default: 4)

=item --base-indent <n>

Sets initial indent level to <n> (default: 1)

=item --dp-to-star-kind <n>

Converts DOUBLE PRECISION to REAL*<dp-to-star-kind> (no default)

=item --dp-to-kind <str>

Converts DOUBLE PRECISION to REAL (KIND=<dp-to-kind>) (no default)

=item --suffix-includes <str>

Adds suffix <str> to file names listed in INCLUDE statements (no default)

=item --verbose

Shows verbose program status information.

=item --version

Shows program version information.

=item --usage

=item --help

Shows standard program information.

=back

=head1 DIAGNOSTICS

Normal termination results in a return value of 0. In case of errors, a
descriptive error message is produced along with a numeric error code.
See EXIT STATUS.

=head1 EXIT STATUS

An exit status code of 0 indicated normal termination. Error codes in
an exit status value are cumulative:

-1) Read failure.
-2) File not found.
-4) Write failure.
-8) Problem during indent.

=head1 CONFIGURATION

None.

=head1 DEPENDENCIES

=over 4

=item *

Carp - Provides saner error messages.

=item *

English - Use human-readable internal variable names.

=item *

Getopt::Long - Process CLI options

=item *

IO::File - Provides cleaner access to files.

=item *

Pod::Usage - Keep usage info consistent with POD documentation

=back

=head1 BUGS AND LIMITATIONS

=over 4

=item *

Some source code confuses the indentation routine.

=item *

The code for handling textual statement labels on DO loops is new
and needs more testing (seems to resolve some issues related to
indenting problem noted above.)

=item *

Strangely-mixed Fortran 77 / Fortran 90 code may confuse f2f.

For simplicity, this code sometimes assumes we are dealing with
a "pure" Fortran 77 input file.  However, many people use mixed
Fortran 77/90/95 code.  This may or may not cause problems.
Give it a try and see; bug reports with sample input which triggers
the bad behavior are cheerfully accepted. :)

=item *

Double precision real variables are considered redundant in
Fortran 90, so this code tries to replace them.  Double precision
variables are "kind=8" on most compilers, but this is not a standard.
Ultimately you will need to determine which 'kind' your compiler thinks
double precision variables are and adjust the code accordingly.
Search for $rh_default->{'dp-to-star-kind'} and $rh_default->{'dp-to-kind'}
in the source code. See below for a reference to MKKIND.F90 which helps
sort out 'kind' issues for various compiler/platform pairs.

=back

=head1 INCOMPATIBILITIES

=over 4

=item * None known.

=back

=head1 TO DO

=over 4

=item *

See TODO comments in source code for potential areas of development

=item *

Isolate and fix remaining indentation problems

=item *

Find and fix bugs

=item *

Allow code to identify F9x free-format code and properly indent it 
(e.g. subsequent runs of this code against a source file should return
the same result; idempotence) 

=back

=head1 AUTHOR

Primary authorship: Colby Lemon, <colbylemon@gmail.com>

Refactoring, documentizing, optionating, and other bastardization courtesy of Bob Apthorpe, <apthorpe+cpan@cynistar.net>

=head1 LICENSE AND COPYRIGHT

Copyright 2003 Colby Lemon <colbylemon@gmail.com>. Some rights reserved. See LICENSE.txt for details.

=head1 SEE ALSO

MKKIND.F90 <http://www.sesp.cse.clrc.ac.uk/Publications/Legacy-Software/Legacy-Software/node38.html>

CONVERT.F90 <http://www.nag.co.uk/nagware/Examples/convert.f90>

=cut
