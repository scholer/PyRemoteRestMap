#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (c) 2015 Rasmus Sorensen, rasmusscholer@gmail.com <scholer.github.io>

##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License

# pylint: disable=C0103,W0142




"""

### Copyright and licenses of the perl library which inspired this python library: ###

Copyright (C) 2003 Peter Blaiklock, pblaiklo@restrictionmapper.org. This module
may be distributed under the same terms as Perl itself.

Source code available from: http://restrictionmapper.org/code.html

   ########


# DESCRIPTION

PyRemoteRestMap.map.Map parses HTML restriction maps generated by sitefind.pl and returns the result as an object.
Sitefind is a CGI program that maps restriction sites.

# Methods

new - Parses the sitefind HTML report and sets attribute values from the result.
Usually called from RemoteRestMap, but can be called by the user with a saved HTML file.

           restriction_map = RemoteRestMap::Map->new($html);

get_next_enz - Returns the next line from the results table as a hash reference.
Returns 0 when done. Take a look at www.restrictionmapper.org for a sample table.
Note that the sorting and filtering options built in to RemoteRestMap/sitefind.pl
can make your life easier. For example, if you find yourself sifting the table for
enzymes with 5' overhangs, just repeat the request with the overhang => ["five_prime"]
parameter to eliminate enzymes with 3' or blunt overhangs.

           while (enz = $restriction_map->get_next_enz) { #returns hash of cut info - enzyme, array of locations, etc.
              enz_name = $enz->{NAME}; # Name of the enzyme
              recognition_site = $enz->{SITE}; # A DNA "regexp" of the enzyme's recognition site
              site_length = $enz->{LENGTH}; # Length of the recognition site in bases
              number_of_cuts = $enz->{CUTNUMBER}; # Number of recognition sites in the sequence
              overhang_type = $enz->{OVERHANG}; # "five_prime", "three_prime" or "blunt"
              cut_locations = @{$enz->{CUTLIST}}; # Reference to an array of cut locations
           }

tab_file - Returns the results table in tab delimited format. Provided as a generic
format for storage, export or whatever.

           file = $restriction_map->tab_file;

noncutters - Returns a list of enzymes that don't cut your sequence. This list is
subject to the original selection criteria that you gave RemoteRestMap (or sitefind).
If an enzyme is not in the table or the noncutters list, it was probably excluded by
these criteria.

           noncutters = $restriction_map->noncutters;

enzyme_list - Returns an alphabetical list of enzymes that cut your sequence.
Note that this may be in a different order than the original table, depending on
the sort options you gave to RemoteRestMap (or sitefind).

           enzymes = $restriction_map->enzyme_list;

cuts - Returns a sorted list of unique positions. Note that this eliminates
duplicate cuts from different enzymes, so that the number of cuts in this list
may be smaller than you expect.

           cuts = $restriction_map->cuts;

total - Returns the total number of rows in the table.

           total = $restriction_map->total;



"""

from .digest import Digest


class MapSites(Digest):
    """
    Represents the result of a "map sites" action on http://www.restrictionmapper.org/
    """
    def __init__(self, ):
        pass

    def standard_headers(self):
        return ["NAME", "SITE", "LENGTH", "CUTNUMBER", "OVERHANG", "CUTLIST"]


    def parse_rows(self, rows, headers):
        if rows is None:
            rows = self.rows
        if headers is None:
            headers = self.headers
        dictrows = (dict(zip(headers, row)) for row in rows)
        for d in dictrows:
            d['CUTPOS'] = [enz.strip() for enz in d[headers[5]].split(",")]
        # Original lib also does an aggregation of cuts and other stuff, but I wait with this until it is asked for.

        return dictrows

    def cuts(self, ):
        """
        Returns a sorted list of unique positions. Note that this eliminates
        duplicate cuts from different enzymes, so that the number of cuts in this list
        may be smaller than you expect.
        """
        return sorted(set(cutpos for row in self.dictrows for cutpos in row['CUTPOS']))
