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



# NAME

RemoteRestMap - an object oriented interface to sitefind.pl


# SYNOPSIS

        from pyremoterestmap import RemoteRestMap
        request = RemoteRestMap(url, dna)
        request.url(new_url)            #get/set method for the URL
        request.sequence(new_sequence)  #get/set method for the DNA sequence

        settings = dict(
           DNAtype       = "circular",
           first         = "site_length",
           second        = "frequency",
           third         = "name",
           enzymetype    = "NEB",
           maxcuts       = "all",
           minlength     = 6,
           isoschizomers = "yes",
           overhang      = ["five_prime", "blunt"],
           enzymelist    = ['BamHI', 'EcoRI']
        )

   restriction_map = request.get_map(\%settings) #returns map object
   while (enz = restriction_map.get_next_enz()) { #returns hash of cut info - enzyme, array of locations, etc.
      enz_name = enz{NAME} # Name of the enzyme
      recognition_site = enz{SITE} # A DNA "regexp" of the enzyme's recognition site
      site_length = enz{LENGTH} # Length of the recognition site in bases
      number_of_cuts = enz{CUTNUMBER} # Number of recognition sites in the sequence
      overhang_type = enz{OVERHANG} # "five_prime", "three_prime" or "blunt"
      cut_locations = @enz{CUTLIST} # Reference to an array of cut locations
   }
   table = restriction_map.tab_file #returns site table in tab delimited format
   noncutters = restriction_map.noncutters #returns array of names of enzymes that don't cut the sequence
   names = restriction_map.enzyme_list #returns alphabetical list of enzyme names
   cuts = restriction_map.cuts #returns unique, ordered list of cut positions

 virtual_digest = request.get_digest({"enzymelist" : enzymes})  #returns digest object.
 while (fragment = virtual_digest.get_next_fragment) {          #returns hash of fragment info
   frag_length = fragment{LENGTH}           # Length of fragment in bases
   five_prime_enz = fragment{FIVEPRIME}     # Name of the enzyme that made the 5' end
   five_prime_pos = fragment{STARTPOS}      # Position of 5' end in original sequence.
   three_prime_enz = fragment{THREEPRIME}   # Name of the enzyme that made the 3' end
   three_prime_pos = fragment{ENDPOS}       # Position of 3' end in original sequence
   dna = fragment{SEQUENCE}                 # DNA sequence of the fragment, 5' to 3'
}
table = virtual_digest.tab_file #returns fragment table in tab delimited format.

# DESCRIPTION

   RemoteRestMap is an object oriented interface to sitefind.pl. Sitefind is a Perl CGI program that maps restriction sites in DNA
   sequences or performs a virtual restriction digest. It provides several options for sorting and filtering output, supports linear
   or circular DNA and can simulate a simultaneous restriction digest. For more information, see www.restrictionmapper.org. RemoteRestMap
   accepts a URL for sitefind, a DNA sequence and a (mostly) optional reference to a hash of parameters. It then validates the parameters,
   supplies default values for missing parameters and submits the request to the URL. The HTML response is parsed and returned to the
   calling program as an object.


## Objects

RemoteRestMap can create 3 objects a request builder, a map object and a digest object. To create a request, pass a URL
and DNA sequence to the RemoteRestMap constructor.

    request = RemoteRestMap(url, dna)

The request object has two attributes url and sequence, both get/set. A good URL for occasional use is
http://www.restrictionmapper.org/cgi-local/sitefind3.pl If you find yourself making frequent heavy use of RemoteRestMap
please install your own copy of sitefind3, it is available under the GPL. More information can be found at www.restrictionmapper.org.
The sequence attribute is just a DNA sequence string. It can contain only numbers, whitespace and the base characters A,C,G and T
(case insensitive).

    request.url(url)
    request.sequence(dna)

RemoteRestMap has two methods, get_map and get_digest. The get_map method returns a map object that is basically a table of data
for each enzyme that cuts the sequence. The get_digest method returns a similar table object for digestion fragments. Both methods
can be passed a hash of parameters to specify sitefind options. These parameters are optional except the enzymelist parameter in get_digest.
For more information see the documentation for these methods and for the RemoteRestMap::Map and RemoteRestMap::Digest modules.

    map = request.get_map(\%parameters)
    digest = request.get_digest({enzymelist : \@enzyme_list})

The map object has six methods, get_next_enz, tab_file, enzyme_list, cuts, cut_number and noncutters.

get_next_enz - returns the next row in the table as a hash or 0 when done. The keys are NAME, SITE, OVERHANG, LENGTH, CUTNUMBER and
CUTLIST. These contain respectively the enzyme name, a DNA "regexp" representing the enzyme's recognition sequence, the type of overhang
the enzyme produces (five_prime, three_prime or blunt), the length of the recognition sequence, the number of times the enzyme cuts the
sequence, and a reference to an array of cut positions.

    while (cuts = restriction_map.get_next_enz) { #returns reference to hash of cut info - enzyme, array of locations, etc.
       enz = cuts.{NAME} # Name of the enzyme
       recognition_site = cuts.{SITE} # A DNA "regexp" of the enzyme's recognition site
       site_length = cuts.{LENGTH} # Length of the recognition site in bases
       number_of_cuts = cuts.{CUTNUMBER} # Number of recognition sites in the sequence
       overhang_type = cuts.{OVERHANG} # five_prime, three_prime or blunt
       cut_locations = @{cuts.{CUTLIST}} # List of cut locations
    }

tab_file - returns the entire enzyme table in tab delimited text format.

    table = restriction_map.tab_file #returns enzyme table in tab delimited format

enzyme_list - returns an alphabetical list of enzymes that cut the sequence. DO NOT use this method if you need to preserve the
original name order, instead, populate your array one entry at a time using get_next_enz.

    enzymes = restriction_map.enzyme_list #returns enzyme list

cuts - returns a unique, ordered list of cut positions. Note that this method will eliminate duplicate cuts produced by different enzymes.

  cuts = restriction_map.cuts #returns list of cut positions

enz_number - returns the number of enzymes that cut your sequence.

  total = restriction_map.total #returns total number of enzymes

non_cutters - returns a list of names of enzymes that do not cut the sequence. Note that these enzymes are subject to the same selection
criteria as the cutters.

  noncutters = restriction_map.noncutters


The get_digest method is similar to get_map. It returns a table object with two methods, get_next_fragment and tab_file. tab_file
returns the fragments table in tab delimited format.

  table = digest.tab_file

get_next_fragment - returns the next row in the fragments table as a hash reference or 0 when finished.
The keys are:
         SEQUENCE - holds the DNA sequence of the fragment.
         LENGTH - the length of the sequence in bases.
         FIVEPRIME - the name of the enzyme that made the 5' cut.
         STARTPOS - the position of the 5' cut in the original sequence.
         THREEPRIME - the name of the enzyme that made the 3' cut.
         ENDPOS - the position of the 3' cut in the original sequence.

    while (fragment = virtual_digest.get_next_fragment) {
        #returns reference to hash of fragment info }
    for fragment in virtual_digest.gen_fragments():
        # process fragments


"""


from .digest import Digest
from .map import Map

import requests
from bs4 import BeautifulSoup


def default_settings(function=None):
    """
    Return default site settings
    """
    defaults = {
                "DNAtype"       : ['linear', 'circular'],
                "first"         : ['frequency', 'overhang', 'name', 'site_length'],
                "second"        : ['overhang', 'frequency', 'name', 'site_length'],
                "third"         : ['name', 'overhang', 'frequency', 'site_length'],
                "enzymetype"    : ['all', 'NEB'],   # all="All Commercial",
                "maxcuts"       : ['all', '0', '1', '2', '3', '4', '5', '10', '20', '30', '40'],
                "minlength"     : ['5', '4', '6', '7', '8'],
                # This settig has been renamed to "Prototypes Only" in the web interface.
                # This is pretty fucked up: "all"->"yes" and "yes"->"no" !!
                "isoschizomers" : ['all', 'yes']
             }
    # Return empty dict to use site defaults:
    if not function:
        return {}
    elif function is 'digest':
        return {
            "DNAtype": "linear",
            "enzymelist": [],   # This is required and must be non-empty when you invoke the request!
            "digest": 1,
            "seqname": "",
        }
    else:
        return {
            "DNAtype": "linear",
            "first": "frequencey", "second": "overhang", "third": "name",   # Sort order
            "enzymetype": "all",
            "maxcuts": "all",
            "minlength": 5,
            "isoschizomers": "no",
            "overhang": ["five_prime", "three_prime", "blunt"],
            "enzymelist": None,        # name = "enzymelist" - has "All Enzymes"
            "digest": 0,
        }

class RemoteRestMap(object):


    def __init__(self, url, sequence, settings=None):
        """
        A good url to use is: http://www.restrictionmapper.org/cgi-local/sitefind3.pl
        """

        self.ValidBases = "ATGC"
        self.AllowNonValidBases = True
        self.StripNonValidBases = True
        self.Url = url
        self.Sequence = sequence
        self.Settings = settings


    @property
    def Sequence(self):
        return self._sequence
    @Sequence.setter
    def Sequence(self, sequence):
        if any(l not in self.ValidBases for l in sequence.upper()):
            if not self.AllowNonValidBases:
                print("Invalid sequence")
                return
        self._sequence = "".join(l for l in sequence.upper()
                                 if l in self.AllowNonValidBases
                                 or not self.StripNonValidBases)

    @property
    def Settings(self):
        if not self._settings:
            return self.default_settings()
        return self._settings
    @Sequence.setter
    def Settings(self, settings):
        if not self._validate_settings(settings):
            print("Invalid settings")
        self._settings = settings

    def get_settings(self, function=None, **kwargs):
        settings = self.default_settings(function)
        settings.update(self.Settings)
        settings.update(kwargs)
        self._validate_settings(settings)   # Should raise an exception if issue.
        return settings

    def _validate_settings(self, settings):
        return True






    def get_digest(self, settings=None):
        """
        Args    : a hash reference of sitefind.pl parameters.

          Possible parameters are:
            DNAtype - Values are "linear" or "circular". This key is optional, default is "linear"
            enzymelist - Value is a reference to an array of enzyme names. This key is required.

        """
        if settings is None:
            settings = self.Settings

        settings = settings.copy()
        settings['digest'] = 1      # Tells sitefind to digest instead of map

        response = self._fetch(self.Url, settings)
        # Probably verify response here...
        digest = Digest(response.text)
        return digest


    def get_map(self, settings):
        """
        Args:
            settings - a dict with parameters for sitefind.pl

        Valid settings/parameters are:
            DNAtype - "linear" or "circular". Specifies DNA conformation.
                Default is "linear".
            first    - "frequency", "name", "overhang" or "site_length".
                Specifies column for first order sorting of output table.
                Default is "frequency".
            second   - "frequency", "name", "overhang" or "site_length".
                Specifies column for second order sorting of output table.
                Default is "overhang".
            third    - "frequency", "name", "overhang" or "site_length".
                Specifies column for third order sorting of output table.
                Default is "name".
            enzymetype - "NEB" or "all".
                Specifies whether all commercial enzymes or only New England Biolabs supplied enzymes
                will be returned.
                Default is "all".
            maxcuts   - "all", "0", "1", "2", "3", "4", "5", "10", "20", "30" or "40".
                Specifies maximum number of cuts per enzyme to return. Enzymes that cut the sequence more
                than this number will not appear in the table. Be sure to use string notation, even for
                numeric values (sorry about that). Use "0" to return noncutters only.
                Default is "all".
            minlength   - 4, 5, 6, 7 or 8.
                Specifies the minimum length of the recognition site.
                Enzymes with sites shorter than this number will not appear in the output. Use numeric notation.
                Default is 5.
            isoschizomers - "yes" or "no".
                Specifies whether to include isoschizomers.
                Default is "no" (prototypes only).
            overhang   -  An array reference containing any combination of "five_prime", "three_prime" and "blunt".
                Specifies which overhangs to include.
                Default is all three.
            enzymelist - A reference to an array of enzyme names. Note that this list overides all other
                         selection criteria (except DNAtype, of course).
                Specifies exactly which enzymes to return.
                Default is no list.
        """
        if settings is None:
            settings = self.Settings
        response = self._fetch(self.Url, settings)
        map_obj = Map(response.text)
        return map_obj


    def _fetch(self, url, settings):
        """
        form action = "cgi-bin/sitefind3.pl
        method = "post"
        onsubmit="return validate_sequence(document.rm_form.sequence.value);
        """
        if url is None:
            url = self.Url
        if settings is None:
            settings = self.Settings

        form = settings.copy()
        user_agent = "PyRemoteRestMap/0.1"
        form['sequence'] = self.Sequence

        # data : is sent in the post request body; params are sent in the query.
        # To debug, use: requests.Request('post', url=url, data=form).prepare().body
        res = requests.post(url, data=form)
        res.raise_for_errors()
        soup = BeautifulSoup(res.text)
        title = soup.find('title').text
        if title == "Error":
            raise ValueError("Error response from %s: %s" % (url, title))

        # time.sleep(3)     # I assume this is in order not to overload the server. Should be done better.

        return res
