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
Copyright and licenes of perl library which inspired this python library:

Copyright (C) 2003 Peter Blaiklock, pblaiklo@restrictionmapper.org. This module
may be distributed under the same terms as Perl itself.

Source code available from: http://restrictionmapper.org/code.html




# DESCRIPTION

PyRemoteRestMap.digest.Digest(html) parses HTML virtual digests generated by sitefind.pl
and returns the result as an object. Sitefind is a CGI program that maps restriction sites.
See http://restrictionmapper.org/

# Methods

get_next_frag - Returns the next line from the results table as a hash reference.
Returns 0 when done. Take a look at www.restrictionmapper.org for a sample table.
See the synopsis for the hash keys and contents.

    frag = digets.get_next_frag()

tab_file - Returns the results table in tab delimited format. Provided as a generic
format for storage, export or whatever.

   file = digest.tab_file()

total - Returns the number of fragments in the digest.

   frag_number = digest.total()



"""




import requests
from bs4 import BeautifulSoup


class Digest(object):
    """
    Represents the result of a "virtual digest".
    NOTE: The fragments produced assumes 100% efficiency of the digest (as implemented on the server side).
    If this is not the case, you may want to extract the restriction sites yourself
    using the map.MapSites data.
    """


    def __init__(self, html, headers=None):

        if headers and headers == "standard":
            self.headers = self.standard_headers()
        else:
            self.headers = headers
        self.root = root = BeautifulSoup(html)
        title = root.find('title')
        if "No Cut Sites" in title.text:
            raise ValueError("No Cut Sites")

        # Table has all the cutters:
        self.rows = None
        self.dictrows = None
        self.Noncutters = None
        self.parse_htmldoc(root)

    def standard_headers(self):
        return ["LENGTH", "START_ENZ", "FIVE_PRIME", "END_ENZ", "THREE_PRIME", "SEQUENCE"]


    def parse_htmldoc(self, root):
        table = root.find('table')
        rows = table.find_all('tr')
        if self.headers is None:
            headerrow = rows.pop(0)
            self.headers = [td.text.strip() for td in headerrow.find_all('td')]
        self.rows = [[td.text.strip() for td in row.find_all('td')] for row in rows]
        # Make dict-list data structure:
        self.dictrows = self.parse_rows()
        # Noncutters:
        self.Noncutters = self.get_noncutters(root)

    def parse_rows(self, rows=None, headers=None):
        if rows is None:
            rows = self.rows
        if headers is None:
            headers = self.headers
        dictrows = [dict(zip(headers, row)) for row in rows]
        return dictrows


    def get_noncutters(self, htmldoc):
        """
        Return a list of non-cutting enzymes from the htmldoc.
        """
        # Wow, this really is the most brain-dead way to deliver data... Give me a JSON interface, ffs.
        bolds = htmldoc.find_all('b')
        noncutters = [elem for elem in bolds if "Noncutters:" in elem.text]
        if noncutters:
            return [enz.strip() for enz in noncutters[0].text.split(',').replace('Noncutters:', '')]



    def gen_fragments(self):
        """
        Neither this, nor get_next_fragment is really needed for python coding, just use:
            for enzdigets

        """
        return (enz for enz in self.dictrows)

    @property
    def total(self):
        return len(self.rows) if self.rows else 0

    def tab_file(self, include_header=True):
        """ Make it easy to make a file from the object. """
        if not self.headers:
            return
        return "\n".join("\t".join(row) for row in
                         (self.headers if include_header and self.headers else []) + self.rows)
