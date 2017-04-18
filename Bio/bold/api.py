# -*- coding: utf-8 -*-
# Copyright 2014 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to invoke the BOLD server over the internet.

This module provides code to work with the BOLD API provided by BOLDSYSTEMS
http://www.boldsystems.org/index.php/Resources

Classes:
Request   Builds the correct URL and performs the HTTP request to the API from
          BOLD.
Response  Parses the returned data from BOLD.

"""
import json
import os
import re
from random import randint
import sys
import warnings
import xml
import xml.etree.ElementTree as ET

from Bio import BiopythonWarning
from Bio import SeqIO
from Bio._py3k import Request as _Request
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import urlencode as _urlencode
from Bio._py3k import _as_string

from . import utils


class Response(object):
    """Accepts and parses results from a call to the BOLD API.

    Parses the data and returns a Response object.

    Attributes:
        items (list or str): Metadata from BOLD after parsing.
        service (str): Alias of the method used to interact with BOLD.

    """
    def _parse_data(self, service, result_string):
        """Parses XML response from BOLD.

        Args:
            service: Alias of the method used to interact with BOLD.
            result_string: XML or JSON string returned from BOLD.

        Returns:
            List of all items as dictionaries.

        """
        self.method = service

        if result_string.strip() == '':
            raise ValueError("BOLD did not return any result.")

        if service == 'call_taxon_search' or service == 'call_taxon_data':
            self._parse_json(result_string)

        if service == 'call_specimen_data' or service == 'call_full_data' or \
                service == 'call_id':
            # Result_string could be data as tab-separated values (tsv)
            # ugly hack for python 2.6 that does not have ET.ParseError
            if sys.version.startswith('2.6'):
                try:
                    self._parse_xml(result_string)
                except xml.parsers.expat.ExpatError:
                    self.items = result_string
            else:
                try:
                    self._parse_xml(result_string)
                except ET.ParseError:
                    self.items = result_string

        if service == 'call_sequence_data':
            self._parse_fasta(result_string)

        if service == 'call_trace_files':
            # file_contents is in binary form
            self.file_contents = result_string

    def _parse_json(self, result_string):
        """Parses JSON response from BOLD.

        Args:
            result_string: JSON string returned from BOLD.

        Returns:
            List of all items as dictionaries.

        Raises:
            ValueError: "BOLD did not return any result."

        """
        items_from_bold = []
        append = items_from_bold.append
        response = json.loads(result_string)
        if hasattr(response, 'items'):
            # Is this a simple JSON and we got only one item?
            simple_json = False
            for i in response.keys():
                res = re.search('^[0-9]+', i)
                if res is None:
                    simple_json = True

            if simple_json is True:
                response = [response]
            for string_id in response:
                item = dict()
                try:
                    json_obj = response[string_id]
                except TypeError:
                    obj = string_id
                    json_obj = obj

                if hasattr(json_obj, 'items'):
                    for k, v in json_obj.items():
                        if k == 'taxid':
                            item['tax_id'] = v
                        elif k == 'taxon':
                            item['taxon'] = v
                        elif k == 'tax_rank':
                            item['tax_rank'] = v
                        elif k == 'tax_division':
                            item['tax_division'] = v
                        elif k == 'parentid':
                            item['parent_id'] = v
                        elif k == 'parentname':
                            item['parent_name'] = v
                        elif k == 'taxonrep':
                            item['taxon_rep'] = v
                        else:
                            item[k] = v
                    append(item)
            self.items = items_from_bold
        else:
            raise ValueError("BOLD did not return any result.")

    def _parse_xml(self, result_string):
        """Parses XML response from BOLD.

        Args:
            result_string: XML string returned from BOLD.

        Returns:
            List of all items as dictionaries.

        """
        items_from_bold = []
        append = items_from_bold.append

        if self.method == 'call_id':
            xml_tag = 'match'
        else:
            xml_tag = 'record'

        root = ET.fromstring(result_string)
        for match in root.findall(xml_tag):
            item = dict()
            fields = [
                # These pairs correspond to convertions of key names from BOLD
                # to friendly versions:
                #
                # (key name from BOLD, friendlier key name)

                # For call_id
                ('ID', 'bold_id'),
                ('sequencedescription', 'sequence_description'),
                ('database', 'database'),
                ('citation', 'citation'),
                ('taxonomicidentification', 'taxonomic_identification'),
                ('similarity', 'similarity'),
                ('specimen/url', 'specimen_url'),
                ('specimen/collectionlocation/country', 'specimen_collection_location_country'),
                ('specimen/collectionlocation/coord/lat', 'specimen_collection_location_latitude'),
                ('specimen/collectionlocation/coord/lon', 'specimen_collection_location_longitude'),

                ('record_id', 'record_id'),
                ('processid', 'process_id'),
                ('bin_uri', 'bin_uri'),
                ('specimen_identifiers/sampleid', 'specimen_identifiers_sample_id'),
                ('specimen_identifiers/catalognum', 'specimen_identifiers_catalog_num'),
                ('specimen_identifiers/fieldnum', 'specimen_identifiers_field_num'),
                ('specimen_identifiers/institution_storing', 'specimen_identifiers_institution_storing'),
                ('taxonomy/identification_provided_by', 'taxonomy_identification_provided_by'),
                ('taxonomy/phylum/taxon/taxID', 'taxonomy_phylum_taxon_id'),
                ('taxonomy/phylum/taxon/name', 'taxonomy_phylum_taxon_name'),
                ('taxonomy/class/taxon/taxID', 'taxonomy_class_taxon_id'),
                ('taxonomy/class/taxon/name', 'taxonomy_class_taxon_name'),
                ('taxonomy/order/taxon/taxID', 'taxonomy_order_taxon_id'),
                ('taxonomy/order/taxon/name', 'taxonomy_order_taxon_name'),
                ('taxonomy/family/taxon/taxID', 'taxonomy_family_taxon_id'),
                ('taxonomy/family/taxon/name', 'taxonomy_family_taxon_name'),
                ('taxonomy/genus/taxon/taxID', 'taxonomy_genus_taxon_id'),
                ('taxonomy/genus/taxon/name', 'taxonomy_genus_taxon_name'),
                ('taxonomy/species/taxon/taxID', 'taxonomy_species_taxon_id'),
                ('taxonomy/species/taxon/name', 'taxonomy_species_taxon_name'),
                ('specimen_details/voucher_type', 'specimen_details_voucher_type'),
                ('specimen_details/voucher_desc', 'specimen_details_voucher_desc'),
                ('specimen_details/extrainfo', 'specimen_details_extra_info'),
                ('specimen_details/lifestage', 'specimen_details_lifestage'),
                ('collection_event/collector', 'collection_event_collector'),
                ('collection_event/collectors', 'collection_event_collectors'),
                ('collection_event/collectiondate', 'collection_event_collection_date'),
                ('collection_event/coordinates/lat', 'collection_event_coordinates_latitude'),
                ('collection_event/coordinates/long', 'collection_event_coordinates_longitude'),
                ('collection_event/exactsite', 'collection_event_exact_site'),
                ('collection_event/country', 'collection_event_country'),
                ('collection_event/province', 'collection_event_province'),
                ('specimen_imagery/media/mediaID', 'specimen_imagery_media_id'),
                ('specimen_imagery/media/caption', 'specimen_imagery_media_caption'),
                ('specimen_imagery/media/metatags', 'specimen_imagery_media_metatags'),
                ('specimen_imagery/media/copyright', 'specimen_imagery_media_copyright'),
                ('specimen_imagery/media/image_file', 'specimen_imagery_media_image_file'),
                ('tracefiles/read/read_id', 'tracefiles_read_read_id'),
                ('tracefiles/read/run_date', 'tracefiles_read_run_date'),
                ('tracefiles/read/sequencing_center', 'tracefiles_read_sequencing_center'),
                ('tracefiles/read/direction', 'tracefiles_read_direction'),
                ('tracefiles/read/seq_primer', 'tracefiles_read_seq_primer'),
                ('tracefiles/read/trace_link', 'tracefiles_read_trace_link'),
                ('tracefiles/read/markercode', 'tracefiles_read_marker_code'),
                ('sequences/sequence/sequenceID', 'sequences_sequence_sequence_id'),
                ('sequences/sequence/markercode', 'sequences_sequence_marker_code'),
                ('sequences/sequence/genbank_accession', 'sequences_sequence_genbank_accession'),
                ('sequences/sequence/nucleotides', 'sequences_sequence_nucleotides'),
            ]
            for field in fields:
                if match.find(field[0]) is not None:
                    key = field[1]
                    matched = match.findall(field[0])
                    if len(matched) == 0:
                        item[key] = None
                    elif len(matched) == 1:
                        item[key] = match.find(field[0]).text
                    elif len(matched) > 1:
                        item[key] = [i.text for i in matched]
            append(item)
        self.items = items_from_bold

    def _parse_fasta(self, result_string):
        """Parses string response from BOLD containing FASTA sequences.

        Args:
            result_string: FASTA sequences as string returned from BOLD.

        Returns:
            List of all items as Biopython SeqRecord objects.

        """
        filename = "tmp_" + str(randint(1, 1000000)) + ".fas"
        with open(filename, "w") as handle:
            handle.write(result_string)
        generator = SeqIO.parse(filename, "fasta")
        self.items = [i for i in generator]
        os.remove(filename)


class Request(object):
    """Constructs a :class:`Request <Request>`. Sends HTTP request.

    Returns:
        A :class:`Response <Response>` object.

    """
    def get(self, service, **kwargs):
        """Does HTTP request to BOLD webservice.

        Args:
            service: The BOLD API alias to interact with.
            kwargs: Paramenters send by users.

        Returns:
            A Response class containing parsed data as attribute `items`.

        """
        params = ''

        if service == 'call_id':
            sequence = utils._prepare_sequence(kwargs['seq'])
            params = _urlencode({'db': kwargs['db'], 'sequence': sequence})

        if service == 'call_taxon_search':
            if kwargs['fuzzy'] is True:
                fuzzy = 'true'
            else:
                fuzzy = 'false'
            params = _urlencode({
                'taxName': kwargs['taxonomic_identification'],
                'fuzzy': fuzzy,
            })

        if service == 'call_taxon_data':
            if kwargs['include_tree'] is False:
                params = _urlencode({
                    'taxId': kwargs['tax_id'],
                    'dataTypes': kwargs['data_type'],
                })
            else:
                params = _urlencode({
                    'taxId': kwargs['tax_id'],
                    'dataTypes': kwargs['data_type'],
                    'includeTree': 'true',
                })

        if service == 'call_specimen_data' or service == 'call_sequence_data' or \
                service == 'call_full_data' or service == 'call_trace_files':
            payload = dict()
            for k, v in kwargs.items():
                if v is not None and k != 'url':
                    payload[k] = v
            params = _urlencode(payload)

        url = kwargs['url'] + "?" + params
        req = _Request(url, headers={'User-Agent': 'BiopythonClient'})
        handle = _urlopen(req)
        response = Response()

        if service == 'call_trace_files':
            binary_result = handle.read()
            response._parse_data(service, binary_result)
        else:
            result = _as_string(handle.read())
            response._parse_data(service, result)
        return response


def request(service, **kwargs):
    """Builds our request based on given arguments. Used internally.

    Args:
        service: The BOLD API alias to interact with. Examples: `call_id`,
                 `call_taxon_search`.
        kwargs: Arguments passed by users when calling our methods.

    Returns:
        Request object with service alias, correct URL and user arguments.

    """
    req = Request()

    if service == 'call_id':
        # User wants the service `call_id`. So we need to use this URL:
        url = "http://boldsystems.org/index.php/Ids_xml"
        return req.get(service=service, url=url, **kwargs)

    if service == 'call_taxon_search':
        url = "http://www.boldsystems.org/index.php/API_Tax/TaxonSearch"
        return req.get(service=service, url=url, **kwargs)

    if service == 'call_taxon_data':
        url = "http://www.boldsystems.org/index.php/API_Tax/TaxonData"
        return req.get(service=service, url=url, **kwargs)

    if service == 'call_trace_files':
        url = "http://www.boldsystems.org/index.php/API_Public/trace"

        args_returning_lots_of_data = ['institutions', 'researchers', 'geo']
        for arg in args_returning_lots_of_data:
            if kwargs[arg] is not None:
                warnings.warn('Requesting ``' + arg + '`` data from BOLD will '
                                                      'possibly return a lot of records and the transfer '
                                                      'of data might take a lot of time to complete as '
                                                      'many Megabytes are expected.',
                              BiopythonWarning
                              )
        return req.get(service=service, url=url, **kwargs)

    if service == 'call_specimen_data':
        url = "http://www.boldsystems.org/index.php/API_Public/specimen"

        args_returning_lots_of_data = ['institutions', 'researchers', 'geo']
        for arg in args_returning_lots_of_data:
            if kwargs[arg] is not None:
                warnings.warn('Requesting ``' + arg + '`` data from BOLD will '
                              'possibly return a lot of records and the transfer '
                              'of data might take a lot of time to complete as '
                              'many Megabytes are expected.',
                              BiopythonWarning
                              )
        return req.get(service=service, url=url, **kwargs)

    if service == 'call_sequence_data':
        url = "http://www.boldsystems.org/index.php/API_Public/sequence"
    elif service == 'call_full_data':
        url = "http://www.boldsystems.org/index.php/API_Public/combined"

    args_returning_lots_of_data = ['institutions', 'researchers', 'geo']
    for arg in args_returning_lots_of_data:
        if kwargs[arg] is not None:
            warnings.warn('Requesting ``' + arg + '`` data from BOLD will '
                                                  'possibly return a lot of records and the transfer '
                                                  'of data might take a lot of time to complete as '
                                                  'many Megabytes are expected.',
                          BiopythonWarning
                          )
    return req.get(service=service, url=url, **kwargs)


def call_id(seq, db):
    """Call the ID Engine API
    http://www.boldsystems.org/index.php/resources/api?type=idengine

    Args:
        seq: DNA sequence string or seq_record object.
        db: The BOLD database of available records. Choices: ``COX1_SPECIES``,'
            ``COX1``, ``COX1_SPECIES_PUBLIC``, ``COX1_L640bp``.

    Returns:
        List of dictionaries containing metadata. One dictionary per BOLD record.

    Examples:

        >>> from Bio import bold
        >>> seq = 'TTTTTGGTATTTGAGCAGGAATAGTAGGAACTTCTCTCAGTTTAATTATTCGAATAGAATTAGGTAATCCAGGTTTCTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCCCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTGTAATTGGAGGATTTGGAAATTGACTAGTTCCCCTAATATTAGGTGCACCTGATATAGCTTTCCCTCGTATAAATAATATAAGATATTGACTACTTCCACCATCTTTAATATTATTAATTTCAAGTAGTATTGTAGAAAATGGAGCTGGAACAGGTTGAACAGTTTACCCCCCTCTTTCCTCTAATATTGCTCATAGAGGAACCTCAGTAGACTTAGCAATTTTTTCTCTTCATTTAGCTGGTATTTCTTCTATTTTAGGAGCTATTAATTTTATTACTACAATTATTAATATACGAGTTAATGGAATATCCTATGATCAAATACCTTTATTTGTTTGAGCTGTTGGAATTACAGCTCTTCTTTTACTTCTTTCTTTACCTGTTTTAGCAGGAGCTATCACAATACTTCTTACAGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGGTGATCCAATTTTATACCAACATTTATTTTGATTTTTTGGTCACCC'
        >>> res = bold.call_id(seq, db='COX1')
        >>> item = res.items[1]
        >>> item['bold_id']  # this is the ID assigned by BOLD
        'GBLN3590-14'

    """
    return request('call_id', seq=seq, db=db)


def call_taxon_search(taxonomic_identification, fuzzy=None):
    """Call the TaxonSearch API
    http://www.boldsystems.org/index.php/resources/api?type=taxonomy#Ideasforwebservices-SequenceParameters

    Args:
        taxonomic_identification: species or any taxon name
        fuzzy: False by default

    Returns:
        List of dictionaries containing metadata. One dictionary per BOLD record.

    Raises:
        ValueError: If `fuzzy` is not True or False.

    Examples:

        >>> from Bio import bold
        >>> taxonomic_identification = 'Euptychia ordinata'
        >>> res = bold.call_taxon_search(taxonomic_identification, fuzzy=False)
        >>> item = res.items[0]  # there can be more than one result
        >>> item['tax_id']
        302603

    """
    if fuzzy is None or fuzzy is False:
        fuzzy = False
    elif fuzzy is True:
        fuzzy = True
    else:
        raise ValueError('Invalid value for ``fuzzy``. Use True or False.')

    return request('call_taxon_search',
                   taxonomic_identification=taxonomic_identification,
                   fuzzy=fuzzy
                   )


def call_taxon_data(tax_id, data_type=None, include_tree=None):
    """Call the TaxonData API. It has several methods to get additional
    metadata.

    Args:
        tax_id: Taxon to get information for.
        data_type: ``basic|all|images``. Default is ``basic``.
        include_tree: Optional. Also returns information for parent taxa. True or
                  False (default).

    Returns:
        List of dictionaries containing metadata for a given taxon.

    Raises:
        ValueError: If `include_tree` is not True or False.

    Examples:

        >>> from Bio import bold
        >>> tax_id = 88899
        >>> res = bold.call_taxon_data(tax_id, data_type='basic,images')
        >>> item = res.items[0]
        >>> item['taxon']
        'Momotus'
        >>> [(i['image'], i['photographer']) for i in item['images']]
        [('BSPBB/MJM_7364_IMG_2240_d+1345758620.JPG', 'Oscar Lopez')]

    """
    if data_type is None:
        # We will use by default data_type='basic'
        data_type = 'basic'

    if include_tree is None or include_tree is False:
        include_tree = False
    elif include_tree is True:
        include_tree = True
    else:
        raise ValueError('Invalid value for ``include_tree``. Use True or False.')

    return request('call_taxon_data', tax_id=tax_id, data_type=data_type,
                   include_tree=include_tree)


def call_specimen_data(taxon=None, ids=None, bin=None, container=None,
                       institutions=None, researchers=None, geo=None,
                       format=None):
    """Call the Specimen Data Retrieval API.

    Args:
        taxon: Taxon name including the ranks: phylum, class, order, family,
               subfamily, genus and species. Example: `taxon='Bos taurus'`.
        ids: Sample ids, process ids, museum ids and field ids. Example:
             `ids='ACRJP618|ACRJP619-11'`.
        bin: BIN stands for Barcode Index number URI. Example: `bin='BOLD:AAA5125'`.
        container: Containers include project codes and dataset codes. Example:
                   `container='DS-EZROM'`.
        institutions: Name of Specimen Storing Sites. Example:
                      `'institutions=Biodiversity Institute of Ontario'`.
        researchers: Collectors and specimen indenfitiers. Example:
                     `researchers='Thibaud Decaens'`.
        geo: Geographic sites such as countries, provinces and states. Example:
             `geo='Alaska'`.
        format: Optional: ``format='tsv'`` will return results a string
                containing data in tab-separated values. If not used, the
                data will be returned as dictionary (default behaviour).

    Raises:
        ValueError: If `format` is not None and not 'tsv'.

    Returns:
        Matching specimen data records as string in TSV format or as list of
        dictionaries.

    Examples:

        >>> from Bio import bold
        >>> bin = 'BOLD:AAE2777'
        >>> res = bold.call_specimen_data(bin=bin)
        >>> class_taxon_names = [item['taxonomy_class_taxon_name'] for item in res.items]
        >>> class_taxon_names[0]
        'Insecta'

    """
    if format is not None and format != 'tsv':
        raise ValueError('Invalid value for ``format``')

    return request('call_specimen_data', taxon=taxon, ids=ids, bin=bin,
                   container=container, institutions=institutions,
                   researchers=researchers, geo=geo, format=format
                   )


def call_sequence_data(taxon=None, ids=None, bin=None, container=None,
                       institutions=None, researchers=None, geo=None,
                       marker=None):
    """Call the Specimen Data Retrieval API.

    Args:
        taxon: Taxon name including the ranks: phylum, class, order, family,
               subfamily, genus and species. Example: `taxon='Bos taurus'`.
        ids: Sample ids, process ids, museum ids and field ids. Example:
             `ids='ACRJP618|ACRJP619-11'`.
        bin: BIN stands for Barcode Index number URI. Example: `bin='BOLD:AAA5125'`.
        container: Containers include project codes and dataset codes. Example:
                   `container='DS-EZROM'`.
        institutions: Name of Specimen Storing Sites. Example:
                      `'institutions=Biodiversity Institute of Ontario'`.
        researchers: Collectors and specimen indenfitiers. Example:
                     `researchers='Thibaud Decaens'`.
        geo: Geographic sites such as countries, provinces and states. Example:
             `geo='Alaska'`.
        marker: Genetic marker code. Example: `marker='COI-5P'`.

    Returns:
        DNA sequences of matching records in FASTA format.

    Examples:

        >>> from Bio import bold
        >>> res = bold.call_sequence_data(taxon='Hermeuptychia', geo='Peru')
        >>> items = res.items
        >>> [item.id for item in items]
        ['GBLN4477-14|Hermeuptychia', 'GBLN4478-14|Hermeuptychia', 'GBLN4479-14|Hermeuptychia']

    """
    return request('call_sequence_data', taxon=taxon, ids=ids, bin=bin,
                   container=container, institutions=institutions,
                   researchers=researchers, geo=geo, marker=marker
                   )


def call_full_data(taxon=None, ids=None, bin=None, container=None,
                   institutions=None, researchers=None, geo=None,
                   marker=None, format=None):
    """Call the Full Data Retrieval API (combined).

    Args:
        taxon: Taxon name including the ranks: phylum, class, order, family,
               subfamily, genus and species. Example: `taxon='Bos taurus'`.
        ids: Sample ids, process ids, museum ids and field ids. Example:
             `ids='ACRJP618|ACRJP619-11'`.
        bin: BIN stands for Barcode Index number URI. Example: `bin='BOLD:AAA5125'`.
        container: Containers include project codes and dataset codes. Example:
                   `container='DS-EZROM'`.
        institutions: Name of Specimen Storing Sites. Example:
                      `'institutions=Biodiversity Institute of Ontario'`.
        researchers: Collectors and specimen indenfitiers. Example:
                     `researchers='Thibaud Decaens'`.
        geo: Geographic sites such as countries, provinces and states. Example:
             `geo='Alaska'`.
        marker: Genetic marker code. Example: `marker='COI-5P'`.
        format: Optional. `format='tsv'`.

    Returns:
        The data is returned as a string in TSV format or list of dicts parsed
        from a XML file.

    Raises:
        ValueError: If `format` is not None or 'tsv'.

    Examples:

        >>> from Bio import bold
        >>> res = bold.call_full_data(taxon='Hermeuptychia', geo='Peru')
        >>> item = res.items[0]
        >>> [item['sequences_sequence_genbank_accession'] for item in res.items]
        ['KF466142', 'KF466143', 'KF466144']

    """
    if format is not None and format != 'tsv':
        raise ValueError('Invalid value for ``format``')

    return request('call_full_data', taxon=taxon, ids=ids, bin=bin,
                   container=container, institutions=institutions,
                   researchers=researchers, geo=geo, marker=marker, format=format
                   )


def call_trace_files(taxon=None, ids=None, bin=None, container=None,
                     institutions=None, researchers=None, geo=None,
                     marker=None):
    """Trace files can be retrieved from BOLD by querying with several parameters.

    Args:
        taxon: Taxon name including the ranks: phylum, class, order, family,
               subfamily, genus and species. Example: `taxon='Bos taurus'`.
        ids: Sample ids, process ids, museum ids and field ids. Example:
             `ids='ACRJP618|ACRJP619-11'`.
        bin: BIN stands for Barcode Index number URI. Example: `bin='BOLD:AAA5125'`.
        container: Containers include project codes and dataset codes. Example:
                   `container='DS-EZROM'`.
        institutions: Name of Specimen Storing Sites. Example:
                      `'institutions=Biodiversity Institute of Ontario'`.
        researchers: Collectors and specimen indenfitiers. Example:
                     `researchers='Thibaud Decaens'`.
        geo: Geographic sites such as countries, provinces and states. Example:
             `geo='Alaska'`.
        marker: Genetic marker code. Example: `marker='COI-5P'`.

    Returns:
        A TAR file consisting of compressed Trace Files (traces in either
        .ab1 or .scf format) along with a file listing the Process ID, taxon and
        marker for each Trace File included.

    Examples:

        >>> from Bio import bold
        >>> res = bold.call_trace_files(taxon='Euptychia mollis',
        ...                             institutions='York University')
        >>> with open("trace_files.tar", "wb") as handle:
        ...     handle.write(res.file_contents)
        4106240

    """
    return request('call_trace_files', taxon=taxon, ids=ids, bin=bin,
                   container=container, institutions=institutions,
                   researchers=researchers, geo=geo, marker=marker
                   )
