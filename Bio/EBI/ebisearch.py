# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provide code to access the EBI Search API.

This module aims to make EBI Search API easy to use. See:
https://www.ebi.ac.uk/ebisearch/overview.ebi

The EBI Search API is a text search engine to access biological data resources
hosted at the European Bioinformatics Institute (EMBL-EBI).

The following code was developed by Berenice Batut (berenice.batut@gmail.com)

References:
Park Y.M., Squizzato S., Buso N., Gur T., Lopez R. (2017)
The EBI search engine: EBI search as a service making biological data accessible
for all.;Nucleic Acids Research, May 2, 2017; DOI: 10.1093/nar/gkx359

"""

import requests

import requires_internet
requires_internet.check()

baseUrl = 'http://www.ebi.ac.uk/ebisearch/ws/rest'
sizeLimit = 100
startLimit = 250000


def get_domain_details(domain):
    """Return a dictionary with the details of a given domain in EBI.

    >>> get_domain_details("allebi")

    domain: domain id in EBI
    """
    url = baseUrl + '/' + domain
    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()


def get_details_and_subdomains(domain, level, verbose=False):
    """Return the details (id and name) for the domain and its subdomains.

    domain: domain id in EBI
    verbose: Boolean to define the printing info
    """
    domain_details = {domain['id']: domain['name']}
    if verbose:
        print("\t" * level + "%s -- %s" % (domain['id'], domain['name']))

    if "subdomains" not in domain:
        return domain_details
    else:
        for subdomain in domain["subdomains"]:
            domain_details.update(
                get_details_and_subdomains(subdomain, level + 1, verbose))
    return domain_details


def get_domains(verbose=False):
    """Return the list of domains in EBI as a dictionary.

    The key being the domain id and the value the domain name.

    verbose: boolean to define the printing info
    """
    allebi = get_domain_details("allebi")
    domain_details = {}
    for domain in allebi["domains"]:
        domain_details.update(get_details_and_subdomains(domain, 0, verbose))
    return domain_details


def print_domain_hierarchy():
    """Print the hierarchy of the domains."""
    get_domains(verbose=True)


def check_domain(domain):
    """Check if a domain exist in EBI.

    domain: id of a domain to check
    """
    if domain not in get_domains(verbose=False):
        err_str = "The domain does not correspond to the id of a known domain"
        err_str += " in EBI. "
        err_str += "The list of EBI domains and their id can be "
        err_str += "accessed with get_domains"
        raise ValueError(err_str)


def get_number_of_results(domain, query):
    """Return the number of results for a query on a specific domain in EBI.

    domain: domain id in EBI
    query: query for EBI
    """
    check_domain(domain)
    url = baseUrl + '/' + domain + '?query=' + query + '&size=0'
    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()['hitCount']


def get_subdomain_fields(domain):
    """Return the fields of a domain and its subdomains.

    domain: domain id in EBI
    """
    fields = {
        "searchable": {},
        "retrievable": {},
        "sortable": {},
        "facet": {},
        "topterms": {}
    }
    if "subdomains" not in domain:
        for field in domain["fieldInfos"]:
            field_desc = field["description"]
            field_id = field["id"]
            for option in field["options"]:
                if option["name"] in fields and option["value"] == "true":
                    fields[option["name"]].setdefault(field_id, field_desc)
    else:
        for subdomain in domain["subdomains"]:
            subdomain_fields = get_subdomain_fields(subdomain)
            fields["searchable"].update(subdomain_fields["searchable"])
            fields["retrievable"].update(subdomain_fields["retrievable"])
            fields["sortable"].update(subdomain_fields["sortable"])
            fields["facet"].update(subdomain_fields["facet"])
            fields["topterms"].update(subdomain_fields["topterms"])
    return fields


def get_fields(domain, verbose=True):
    """Return the fields (for different type) of a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    check_domain(domain)
    domain_details = get_domain_details(domain)
    fields = {
        "searchable": {},
        "retrievable": {},
        "sortable": {},
        "facet": {},
        "topterms": {}
    }
    for domain in domain_details["domains"]:
        subdomain_fields = get_subdomain_fields(domain)
        fields["searchable"].update(subdomain_fields["searchable"])
        fields["retrievable"].update(subdomain_fields["retrievable"])
        fields["sortable"].update(subdomain_fields["sortable"])
        fields["facet"].update(subdomain_fields["facet"])
        fields["topterms"].update(subdomain_fields["topterms"])

    if verbose:
        for field_type in fields:
            print("%s" % (field_type))
            for field in fields[field_type]:
                print("\t%s" % (field))
    return fields


def get_specific_fields(domain, field_type, verbose=True):
    """Return the fields of a given type for a specific domain in EBI.

    domain: domain id in EBI
    field_type: type of field to extract (searchable, retrievable, sortable,
    facet, topterms)
    verbose: boolean to define the printing info
    """
    field_types = [
        "searchable", "retrievable", "sortable", "facet", "topterms"]
    if field_type not in field_types:
        err_str = "The type of field to extract must be 'searchable', "
        err_str += "'retrievable', 'sortable', 'facet' or 'topterms'"
        raise ValueError(err_str)

    fields = get_fields(domain, verbose=False)[field_type]
    if verbose:
        print("%s fields for %s" % (field_type, domain))
        for field in fields:
            print("%s" % (field))
    return fields


def get_searchable_fields(domain, verbose=True):
    """Return the searchable fields for a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    return get_specific_fields(domain, "searchable", verbose)


def get_retrievable_fields(domain, verbose=True):
    """Return the retrievable fields of a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    return get_specific_fields(domain, "retrievable", verbose)


def get_sortable_fields(domain, verbose=True):
    """Return the sortable fields of a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    return get_specific_fields(domain, "sortable", verbose)


def get_facet_fields(domain, verbose=True):
    """Return the facet fields of a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    return get_specific_fields(domain, "facet", verbose)


def get_topterms_fields(domain, verbose=True):
    """Return the topterms fields of a specific domain in EBI.

    domain: domain id in EBI
    verbose: boolean to define the printing info
    """
    return get_specific_fields(domain, "topterms", verbose)


def check_order(order):
    """Check if an order is either ascending or descending.

    order: order to check
    """
    if order != "ascending" and order != "descending":
        err_str = "Order value must be either 'ascending' or 'descending'"
        raise ValueError(err_str)


def check_retrievable_fields(fields, domain):
    """Check if the fields are retrievable for a domain.

    field: field to check
    domain: domain id in EBI (accessible with get_domains)
    """
    retrievable_fields = get_retrievable_fields(domain, verbose=False)
    for field in fields.split(","):
        if field not in retrievable_fields:
            err_str = "The field %s does not correspond " % (field)
            err_str += "to a retrievable field for the domain. "
            err_str += "The list of retrievable fields for a domain can be "
            err_str += "accessed with get_retrievable_fields"
            raise ValueError(err_str)


def check_sortable_field(field, domain):
    """Check if a field is sortable for the domain.

    field: field to check
    domain: domain id in EBI (accessible with get_domains)
    """
    sortable_fields = get_sortable_fields(domain, verbose=False)
    if field not in sortable_fields:
        err_str = "The field %s is not a sortable field " % (field)
        err_str += "for the domain %s" % (domain)
        err_str += "The list of sortable fields for the domain can be "
        err_str += "accessed with get_sortable_fields"
        raise ValueError(err_str)


def check_facet_field(field, domain):
    """Check if a field is facet for the domain.

    field: field to check
    domain: domain id in EBI (accessible with get_domains)
    """
    facet_fields = get_facet_fields(domain, verbose=False)
    if field not in facet_fields:
        err_str = "The field %s is not a facet field " % (field)
        err_str += "for the domain %s" % (domain)
        err_str += "The list of sortable fields for the domain can be "
        err_str += "accessed with get_facet_fields"
        raise ValueError(err_str)


def check_topterms_field(field, domain):
    """Check if a field is topterm for the domain.

    field: field to check
    domain: domain id in EBI (accessible with get_domains)
    """
    topterms_fields = get_topterms_fields(domain, verbose=False)
    if field not in topterms_fields:
        err_str = "The field %s is not a topterms field " % (field)
        err_str += "for the domain %s" % (domain)
        err_str += "The list of sortable fields for the domain can be "
        err_str += "accessed with get_topterms_fields"
        raise ValueError(err_str)


def check_size(size):
    """Check that the size is lower than a given limit.

    size: value to check
    limit: threshold
    """
    if size > sizeLimit:
        err_str = "Size (number of entries to retrieve) must be lower "
        err_str += "than %s" % (sizeLimit)
        raise ValueError(err_str)


def check_start(start):
    """Check that the start is lower than a given limit.

    start: value to check
    limit: threshold
    """
    if start > startLimit:
        err_str = "Start (index of the first entry in the results) "
        err_str += "must be lower than %s" % (startLimit)
        raise ValueError(err_str)


def get_domain_search_results(
    domain, query, fields, size=15, start=0, order=None, sortfield=None,
    sort=None, fieldurl=False, viewurl=False, facets=None, facetfields=None,
    facetcount=None, facetsdepth=None
):
    """Return the results for a query on a specific domain in EBI.

    domain: domain id in EBI (accessible with get_domains)
    query: query for EBI (the searchable fields can be accessed with
    get_searchable_fields)
    fields: fields to retrieve for the query (the fields can be accessed with
    get_retrievable_fields), separated by comma
    size: number of entries to retrieve
    start: index of the first entry in the results
    order: order to sort the results (ascending or descending), should come
    along with "sortfield" and not allowed to use with "sort" parameters
    sortfield: single field identifier to sort on (the fields can be accessed
    with get_sortable_fields)
    sort: comma separated values of sorting criteria with field_id:order
    (e.g. boost:descending,length:descending), should not be used in
    conjunction with any of 'sortfield' and 'order' parameters
    fieldurl: boolean to indicate whether field links are included (the
    returned links mean direct URLs to the data entries in original portals)
    viewurl: boolean to indicate whether other view links (than fieldurl) on an
    entry are included
    facets: comma separated values of facet selections to apply on search
    results with facet_id:facet_value (e.g. keywords:Glycolysis). The facet id
    can be accessed with get_facet_fields
    facetfields: comma separated values of field identifiers associated with
    facets to retrieve
    facetcount: number of facet values to retrieve
    facetsdepth: number of level in the hierarchy to retrieve
    """
    url = baseUrl + '/'

    check_domain(domain)

    url += domain
    url += '?query=' + query

    check_retrievable_fields(fields, domain)
    url += '&fields=' + fields

    result_nb = get_number_of_results(domain, query)
    check_size(size)
    if size > result_nb:
        err_str = "Size (number of entries to retrieve) must be lower "
        err_str += "than the number of expected results for the query"
        raise ValueError(err_str)
    url += '&size=%s' % (size)

    check_start(start)
    if start > result_nb:
        err_str = "Start (index of the first entry in the results) must "
        err_str += "be lower than the number of expected results for the "
        err_str += "query"
        raise ValueError(err_str)
    url += '&start=%s' % (start)

    if order is not None or sortfield is not None:
        if sortfield is None:
            err_str = "Order should come along with 'sortfield' parameters"
            raise ValueError(err_str)
        if order is None:
            err_str = "Sortfield should come along with 'order' parameters"
            raise ValueError(err_str)
        check_order(order)
        check_sortable_field(sortfield, domain)
        url += '&order=%s' % (order)
        url += '&sortfield=%s' % (sortfield)

    if sort is not None:
        if order is not None or sortfield is not None:
            err_str = "Sort should not come along with 'sortfield' or 'order' "
            err_str += "parameters"
            raise ValueError(err_str)
        for criteria in sort.split(","):
            split_criteria = criteria.split(":")
            check_sortable_field(split_criteria[0], domain)
            check_order(split_criteria[1])

    if fieldurl:
        url += '&fieldurl=true'
    else:
        url += '&fieldurl=false'

    if viewurl:
        url += '&viewurl=true'
    else:
        url += '&viewurl=false'

    if facets is not None:
        for facet in facets.split(","):
            check_facet_field(facet.split(":")[0], domain)
        url += '&facets=%s' % (facets)

    if facetfields is not None:
        url += '&facetfields=%s' % (facetfields)

    if facetcount is not None:
        url += '&facetcount=%s' % (facetcount)

    if facetsdepth is not None:
        url += '&facetsdepth=%s' % (facetsdepth)

    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()['entries']


def get_all_domain_search_results(
    domain, query, fields, order=None, sortfield=None, sort=None,
    fieldurl=False, viewurl=False, facets=None, facetfields=None,
    facetcount=None, facetsdepth=None
):
    """Return the all the results for a query on a specific domain in EBI.

    domain: domain id in EBI (accessible with get_domains)
    query: query for EBI (the searchable fields can be accessed with
    get_searchable_fields)
    fields: fields to retrieve for the query (the fields can be accessed with
    get_retrievable_fields), separated by comma
    order: order to sort the results (ascending or descending), should come
    along with "sortfield" and not allowed to use with "sort" parameters
    sortfield: single field identifier to sort on (the fields can be accessed
    with get_sortable_fields)
    sort: comma separated values of sorting criteria with field_id:order
    (e.g. boost:descending,length:descending), should not be used in
    conjunction with any of 'sortfield' and 'order' parameters
    fieldurl: boolean to indicate whether field links are included (the
    returned links mean direct URLs to the data entries in original portals)
    viewurl: boolean to indicate whether other view links (than fieldurl) on an
    entry are included
    facets: comma separated values of facet selections to apply on search
    results with facet_id:facet_value (e.g. keywords:Glycolysis). The facet id
    can be accessed with get_facet_fields
    facetfields: comma separated values of field identifiers associated with
    facets to retrieve
    facetcount: number of facet values to retrieve
    facetsdepth: number of level in the hierarchy to retrieve
    """
    result_nb = get_number_of_results(domain, query)
    quotient = int(result_nb / float(sizeLimit))
    start = 0
    all_results = []
    for i in range(quotient):
        start = sizeLimit * i
        all_results.extend(get_domain_search_results(
            domain=domain,
            query=query,
            fields=fields,
            size=sizeLimit,
            start=start,
            order=order,
            sortfield=sortfield,
            sort=sort,
            fieldurl=fieldurl,
            viewurl=viewurl,
            facets=facets,
            facetfields=facetfields,
            facetcount=facetcount,
            facetsdepth=facetsdepth))
    if (result_nb % 100) > 0:
        start = sizeLimit * quotient
        remainder = result_nb - start
        all_results.extend(get_domain_search_results(
            domain=domain,
            query=query,
            fields=fields,
            size=remainder,
            start=start,
            order=order,
            sortfield=sortfield,
            sort=sort,
            fieldurl=fieldurl,
            viewurl=viewurl,
            facets=facets,
            facetfields=facetfields,
            facetcount=facetcount,
            facetsdepth=facetsdepth))
    return all_results


def get_entries(domain, entryids, fields, fieldurl=False, viewurl=False):
    """Return content of entries on a specific domain in EBI.

    domain: domain id in EBI (accessible with get_domains)
    entryids: comma seperated values of entry identifiers
    fields: fields to retrieve for the query (the fields can be accessed with
    get_retrievable_fields), separated by comma
    fieldurl: boolean to indicate whether field links are included (the
    returned links mean direct URLs to the data entries in original portals)
    viewurl: boolean to indicate whether other view links (than fieldurl) on an
    entry are included
    """
    url = baseUrl + '/'

    check_domain(domain)
    url += domain

    url += '/entry/' + entryids

    check_retrievable_fields(fields, domain)
    url += '?fields=' + fields

    if fieldurl:
        url += '&fieldurl=true'
    else:
        url += '&fieldurl=false'

    if viewurl:
        url += '&viewurl=true'
    else:
        url += '&viewurl=false'

    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()['entries']


def get_field_topterms(
    domain, fieldid, size=15, excludes=None, excludesets=None
):
    """Return a list of top terms in a field of a specific domain in EBI.

    domain: domain id in EBI (accessible with get_domains)
    fieldid: field id (the fields can be accessed with get_topterms_fields)
    size: number of entries to retieve
    excludes: comma separated values of terms to be excluded
    excludesets: comma separated values of stop-word sets to be excluded
    """
    url = baseUrl + '/'

    check_domain(domain)
    url += domain

    check_topterms_field(fieldid, domain)
    url += '?fieldid=' + fieldid

    check_size(size)
    url += '&size=%s' % (size)

    if excludes is not None:
        url += '&excludes=%s' % (excludes)

    if excludesets is not None:
        url += '&excludesets=%s' % (excludesets)


def get_number_of_morelikethis(domain, entryid):
    """"Return the number of entries similar to an entry of a domain in EBI.

    domain: domain id in EBI
    entryid: entry id
    """
    url = baseUrl + '/'

    check_domain(domain)
    url += domain
    url += '/entry/' + entryid
    url += '/morelikethis'
    url += '?size=0'
    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()['hitCount']


def get_morelikethis(
    domain, entryid, size=15, start=0, fields=None, fieldurl=False,
    viewurl=False, mltfields=None, mintermfreq=None, mindocfreq=None,
    maxqueryterm=None, excludes=None, excludesets=None
):
    """Return a list of similar entries to an entry of a specific domain in EBI.

    domain: domain id in EBI (accessible with get_domains)
    entryid: entry id
    size: number of entries to retieve
    start: index of the first entry in the results
    fields: field id (the fields can be accessed with get_topterms_fields)
    fieldurl: boolean to indicate whether field links are included (the
    returned links mean direct URLs to the data entries in original portals)
    viewurl: boolean to indicate whether other view links (than fieldurl) on an
    entry are included
    mltfields: comma separated values of field identifiers to be used for
    generating a query
    mintermfreq: minimum term frequency (any terms whose frequency is below
    this value will be ignore from a base entry)
    mindocfreq: maximum document frequency (any terms which occur in at least
    this number of entries will be ignored)
    maxqueryterm: maximum number of query terms that will be included in any
    generated query (max. 25)
    excludes: comma separated values of terms to be excluded
    excludesets: comma separated values of stop-word sets to be excluded
    """
    url = baseUrl + '/'

    check_domain(domain)
    url += domain

    url += '/entry/' + entryid
    url += '/morelikethis'

    result_nb = get_number_of_morelikethis(domain, entryid)
    check_size(size)
    if size > result_nb:
        err_str = "Size (number of entries to retrieve) must be lower "
        err_str += "than the number of expected results for the query"
        raise ValueError(err_str)
    url += '?size=%s' % (size)

    check_start(start)
    if start > result_nb:
        err_str = "Start (index of the first entry in the results) must "
        err_str += "be lower than the number of expected results for the "
        err_str += "query"
        raise ValueError(err_str)
    url += '&start=%s' % (start)

    if fields is not None:
        url += '&fields=%s' % (fields)

    if fieldurl:
        url += '&fieldurl=true'
    else:
        url += '&fieldurl=false'

    if viewurl:
        url += '&viewurl=true'
    else:
        url += '&viewurl=false'

    if mltfields is not None:
        url += '&fields=%s' % (mltfields)
    if mintermfreq is not None:
        url += '&fields=%s' % (mintermfreq)
    if mindocfreq is not None:
        url += '&fields=%s' % (mindocfreq)
    if maxqueryterm is not None:
        url += '&fields=%s' % (maxqueryterm)
    if excludes is not None:
        url += '&fields=%s' % (excludes)
    if excludesets is not None:
        url += '&fields=%s' % (excludesets)

    r = requests.get(
        url,
        headers={"accept": "application/json"})
    r.raise_for_status()
    return r.json()['entries']
