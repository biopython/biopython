"""Parse a time or date string.

Converts a strftime/strptime-like format string into a Martel regular
expression (either as a string or an Expression).

Example use:
  >>> from Martel import Time
  >>> from xml.sax import saxutils
  >>> format = Time.make_expression("%(Jan)-%(day)-%(YYYY)\n")
  >>> parser = format.make_parser()
  >>> parser.setContentHandler(saxutils.XMLGenerator())
  >>> parser.parseString("OCT-31-2021\\n")
  <?xml version="1.0" encoding="iso-8859-1"?>
  <month type="short">OCT</month>-<day type="numeric">31</day>-<year type="long">2021</year>
  >>>
  

Times and dates come up often in parsing.  It's usually pretty easily
to write a syntax for them on the fly.  For example, suppose you want
to parse a date in the format
  YYYY-MM-DD

as in "1985-03-26".  One pattern for that is
  "\\d{4}-\\d{2}-\\d{2}"

To get the individual fields in Martel requires group names.
  "(?P<year>\\d{4})-(?P<month>\\d{2})-(?P<day>\\d{2})"

If you want some minimal verification (eg, to help make sure you
haven't accidentally swapped the day and month fields) you need to
tighten down on what values are allowed, as in
  "(?P<year>\\d{4})-(?P<month>0[1-9]|1[012])-(?P<day>0[1-9]|[12][0-9]|3[01])"

The more you write, the more the likelihood of making a mistake, the
more the chance different format definitions use different patterns,
the harder it is to understand what's going on.

This module helps by providing a set of standard definitions for the
different terms needed in parsing dates, and a way to generate those
definitions from a relatively easy to understand format string.

The syntax of the format string is based on that used by the standard
unix strftime/strptime functions, with terms taken from the POSIX and
GNU documentation plus some experimentation.  These terms are in the
form "%c" where "c" is a single character.  It's hard to remember
everything through a not always mnemonic single character code, so
Martel.Time adds a new syntax of the form "%(word)" where word can be
one of the single characters, or a multicharacter word.  For example,
"%(Mon)" is identical to "%a" but easier to understand.

The complete list of definitions is given below.

The lowest-level terms (like "year", but excluding terms like "%D"
which expand to other terms) are inside of named groups, which
generate the element tag and attributes when used for Martel.

For example, "%m" generates the pattern used for a month, from "01" to
"12".  The name of the group is "month" and it has a single attribute
named "type" with value "numeric".  (All "numeric" types can be parsed
with Python's 'int' function.)  The default pattern made from "%m" is

    (?P<month?type=numeric>(0[1-9]|1[012]))

and when parsed against a month value, like "05", produces

    <month type="numeric">05</month>

The "type" attribute is used because values which mean the same thing
can be represented in different formats.  The month "January" can be
represented with the word "January" (type = "long"), "Jan" (type =
"short"), "01" (type = "numeric"), "1" (type = "numeric"), or " 1"
(type = "numeric").  [Note: It is possible that subtypes may be added
in the future to distinguish between these different numeric cases.]


FUNCTIONS:

There are two public functions -- "make_pattern" and
"make_expression".  Both take the same parameters and return a regular
expression.

  make_pattern(format, tag_format = "%s") -- returns the expression
       as a pattern string

  make_expression(format, tag_format = "%s") -- returns the expression
       as a Martel.Expression data structure (which can be used to
       make a parser)

The first parameter, "format", is the time format string already
discussed.  Some examples are:

  >>> from Martel import Time
  >>> Time.make_pattern("%y")
  '(?P<year?type=short>\\\\d{2})'
  >>> Time.make_pattern("%H:%M")
  '(?P<hour?type=24-hour>([01][0-9]|2[0-3]))\\\\:(?P<minute?type=numeric>[0-5][0-9])'
  >>>

The second parameter is used if you want to change the tag name.  For
example, instead of "year" you may want "year-modified" or
"start-year" -- or you may not want a tag at all.

For each term, the tag name ("year", "month", etc.) is %'ed with the
tag_format string.  The default string is "%s" which effectively says
to keep the name unchanged.  Here are a couple examples which use a
different string.

  >>> from Martel import Time
  >>> Time.make_pattern("%(year)", "%s-modified")
  '(?P<year-modified?type=any>([0-9]{2}([0-9]{2})?))'
  >>> Time.make_pattern("%(year)", "start-%s")
  '(?P<start-year?type=any>([0-9]{2}([0-9]{2})?))'
  >>> Time.make_pattern("%(year)", None)
  '([0-9]{2}([0-9]{2})?)'
  >>>
  
The tag_format is used for every tag name, which lets you modify
several values at once.  You can even pass in an object which
implements the __mod__ method to make more drastic changes to the
name.

  >>>  Time.make_pattern("%H:%M", "%s-created")
  '(?P<hour-created?type=24-hour>([01][0-9]|2[0-3]))\\\\:(?P<minute-created?type=numeric>[0-5][0-9])'
  >>> class Upcase:
  ...     def __mod__(self, name):
  ...         return name.upper()
  ...
  >>> Time.make_pattern("%H:%M", Upcase())
  '(?P<HOUR?type=24-hour>([01][0-9]|2[0-3]))\\:(?P<MINUTE?type=numeric>[0-5][0-9])'
  >>>

BUGS:
Only the "C" locale (essentialy, US/English) is supported.  Field
widths (as in "%-5d") are not supported.

There is no way to change the element attributes.  I'm not sure this
is a bug.

         ====  Table of Date/Time Specifiers ====

%a is replaced by the pattern for the abbreviated weekday name.
  Pattern: (Mon|Tue|Wed|Thu|Fri|Sat|Sun)
     (the real pattern is case insensitive)
  Example: "Wed" "FRI"
  Element name: "weekday"
  Element attributes: "type" = "short"
  Note: %(Mon) is the same as %a
  
%A is replaced by the pattern for the full weekday name.
  Pattern: (Monday|Tuesday|Wednesday|Thursday|Friday|Saturday|Sunday)
     (the real pattern is case insensitive)
  Example: "Thursday" "SUNDAY"
  Element name: "weekday"
  Element attributes: "type" = "long"
  Note: %(Monday) is the same as %a

%b is replaced by the the pattern for the abbreviated month name.
  Pattern: (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)
     (the real pattern is case insensitive)
  Example: "Oct" "AUG"
  Element name: "month"
  Element attributes: "type" = "short"
  Note: %(Jan) is the same as %b
  
%B is replaced by the pattern for the full month name.
  Pattern: (January|February|March|April|May|June|July|August|
            September|October|November|December)
     (the real pattern is case insensitive)
  Example: "August", "MAY"
  Element name: "month"
  Element attributes: "type" = "long"
  Note: %(January) is the same as %B

%c is replaced by the pattern for the US 24-hour date and time
   representation.
  Pattern: same as "%a %b %e %T %Y"
  Example: "Wed Dec 12 19:57:22 2001"
  Element: only uses names and attributes of the individual terms

%C is replaced by the pattern for the century number (the year divided
   by 100 and truncated to an integer) as a decimal number [00-99].
  Pattern: "[0-9][0-9]"
  Example: "19" for the years 1900 to 1999
  Element name: "century"
  Element attributes: "type" = "numeric"

%d is replaced by the pattern for a day of the month as a decimal
   number [01,31].
  Pattern: (0[1-9]|[12][0-9]|3[01])
  Example: "01", "12"
  Element name: "day"
  Element attributes: "type": "numeric"
  Note: "%d" does not include " 1" or "1".  If you also want to allow
    those then use "%(day)"
  
%D same as the pattern for "%m/%d/%y".
  Pattern: see "%m/%d/%y".
  Example: "12/13/01"
  Element: only uses names and attributes of the individual terms
  
%e is replaced by the pattern for a day of the month as a decimal
   number [1,31]; a single digit is preceded by a space.
  Pattern: "( [1-9]|[12][0-9]|3[01])"
  Example: " 1", "31"
  Element name: "day"
  Element attributes: "type" = "numeric"
  Note: "%e" does not include "01" or "1".  If you also want to allow
    those then use "%(day)"
  
%F same as the pattern for "%Y-%m-%d".
  Pattern: see "%Y-%m-%d".
  Example: "2001-12-21"
  Element: only uses names and attributes of the individual terms

%g ISO 8601 2-digit (like %G but without the century) (00-99)
  Pattern: [0-9][0-9]
  Example: "00"
  Element name: "century"
  Element attributes: "type" = "ISO8601"

%G Pattern for the ISO 8601 year with century as a decimal number.
   The 4-digit year corresponding to the ISO week number (see %V).
   This has the same format and value as %y, except that if the ISO
   week number belongs to the previous or next year, that year is
   used instead. (TZ)
  Pattern: [0-9][0-9][0-9][0-9]
  Example: "1954" "2001"
  Element name: "year"
  Element attributes: "type" = "ISO8601"

%h (DEPRECATED) same as %b.  
  Pattern: (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)
  Example: "Feb"
  Element name: "month"
  Element attributes: "type" = "short"
  Note: %(Jan) is the same as %b is the same as %h
  
%H is replaced by the pattern for the hour on a 24-hour clock, as
   a decimal number [00,23].
  Pattern: ([01][0-9]|2[0-3])
  Example: "00", "01", "23"
  Element name: "hour"
  Element attributes: "type" = "24-hour"
  Note: This does not allow single digit hours like "1".  If you also
    want to include those, use %(24-hour)
  
%I is replaced by the pattern for the hour on a 12-hour clock, as
   a decimal number [01,12].
  Pattern: (0[0-9]|1[012])
  Example: "01", "12"
  Element name: "hour"
  Element attributes: "type" = "12-hour"
  Note: This does not allow single digit hours like "1".  If you also
    want to include those, use %(12-hour)
  
%j is replaced by the pattern for day of the year as a decimal
   number.  First day is numbered "001" [001,366].
  Pattern: "([12][0-9][0-9]|3([012345][0-9]|6[0-6])|0(0[1-9]|[1-9][0-9]))"
  Example: "001", "092", "362"
  Element name: "year_day"
  Element attributes: "type": "1"

%k is replaced by the pattern for the hour on a 24-hour clock, as a
   decimal number (range 0 to 23); single digits are preceded by a
   blank.
  Pattern: "( [0-9]|1[0-9]|2[0123])"
  Example: " 1", "10", "23"
  Element name: "hour"
  Element attributes: "type" = "24-hour"
  Note: This does not allow single digit hours like "1" or hours which
    start with an "0" like "03".  If you also want to include those,
    use %(24-hour).  See also %H.

%l is replaced by the pattern for the hour on a 12-hour clock, as a
   decimal number (range 1 to 12); single digits are preceded by a
   blank.
  Pattern: "( [0-9]|1[012])"
  Example: " 1", "10", "12"
  Element name: "hour"
  Element attributes: "type" = "12-hour"
  Note: This does not allow single digit hours like "1" or hours which
    start with an "0" like "03".  If you also want to include those,
    use %(12-hour).  See also %I.

%m is replaced by the pattern for the month as a decimal number [01,12].
  Pattern: "(0[1-9]|1[012])"
  Example: "01", "09", "12"
  Element name: "month"
  Element attributes: "type" = "numeric"
  Note: This does not allow single digit months like "1" or months which
    start with an space like " 3".  If you also want to include those,
    use %(month).  See also %(DD), which is an alias for %m.
  
%M is replaced by the pattern for the minute as a decimal number [00,59].
  Pattern: "[0-5][0-9]"
  Example: "00", "38"
  Element name: "minute"
  Element attributes: "type" = "numeric"
  Note: this is the same as %(minute)
  
%n is replaced by the pattern for the newline character.
  Pattern: "\\n"
  Note: you shouldn't need to use this

%p is replaced by the case insensitive pattern for "AM" or "PM"
  Pattern: "([AaPp][Mm])"
  Example: "AM", "pm"
  Element name: "ampm"
  Element attributes: no attributes
  Note: this doesn't allow "a.m." or "P.M."
  
%P is identical to "%p" (they have slightly different meanings for output)
  Pattern: "([AaPp][Mm])"
  Example: "am", "PM"
  Element name: "ampm"
  Element attributes: no attributes
  Note: this doesn't allow "a.m." or "P.M."

%r is equivalent to "%I:%M:%S %p".
  Pattern: see the patterns for the individual terms
  Example: "07:57:22 PM"
  Element: only uses names and attributes of the individual terms
  
%R is the pattern for the 24 hour notation "%H:%M".
  Pattern: see the patterns for the individual terms
  Example: "19:57"
  Element: only uses names and attributes of the individual terms
  
%s is pattern for a decimal number of seconds (Unix timestamp)
  Pattern: "[0-9]+"
  Example: "1008205042"
  Element name: "timestamp"
  Element attributes: no attributes
  
%S is replaced by the pattern for the second as a decimal number
   Can take values from "00" to "61" (includes double leap seconds).
  Pattern: "([0-5][0-9]|6[01])"
  Example: "03", "25"
  Element name: "second"
  Element attributes: "type" = "numeric"
  Note: This is the same as %(second)
  
%t is replaced by a tab character. (plat-spec)
  Pattern: "\\t"
  Note: You shouldn't need to use this.
  
%T is identical to the 24-hour time format "%H:%M:%S".
  Pattern: see the patterns for the individual terms
  Example: "19:57:22"
  Element: only uses names and attributes of the individual terms
  
%u is replaced by the pattern for the weekday as a decimal number
   [1,7], with "1" representing Monday.
  Pattern: "[1-7]"
  Example: "4" (which is Thursday)
  Element name: "weekday"
  Element attributes: "type" = "Monday1"
  Note: See also %w, which has a type of "Sunday0"
  
%U is replaced by the pattern for the week number of the year (Sunday
   as the first day of the week) as a decimal number [00,53].  In
   other words, this is the number of Sundays seen so far in the year.
  Pattern: "([0-4][0-9]|5[0-3])"
  Example: "04", "26"
  Element name: "week_number"
  Element attributes: "type" = "Sunday_count"
  Note: See also %V and %W
  
%V is replaced by the pattern for the week number of the year (Monday
   as the first day of the week) as a decimal number [01,53]. This is
   used for week numbers where if the week containing 1 January has four
   or more days in the new year, then it is considered week 1. (Otherwise,
   it is the last week of the previous year, and the next week is week 1.)
  Pattern: "(0[1-9]|[1-4][0-9]|5[0-3])"
  Example: "04", "33"
  Element name: "week_number"
  Element attributes: "type" = "type_V" (Got a better short name?)
  Note: See also %U and %W.  I don't know when to use this.
  
%w is replaced by pattern for the the weekday as a decimal number [0,6],
   with 0 representing Sunday.
  Pattern: "[0-6]"
  Example: "6"
  Element name: "weekday"
  Element attributes: "type" = "Sunday0"
  Note: See also %u, which has a type of "Monday1"
  
%W is replaced by the pattern for the week number of the year (Monday
   as the first day of the week) as a decimal number [00,53]. All days
   in a new year preceding the first Monday are considered to be in
   week 0.  In other words, this is the number of Mondays seen so far
   in the year.
  Pattern: "([0-4][0-9]|5[0-3])"
  Example: "00", "49"
  Element name: "week_number"
  Element attributes: "type" = "Monday_count"
  Note: See also %U and %V.

%x is the same as "%D", which is "%m/%d/%y".
  Pattern: see the patterns for the individual terms
  Example: "12/13/99"
  Element: only uses names and attributes of the individual terms
  
%X is the same as "%T", which is "%H:%M:%S".
  Pattern: see the patterns for the individual terms
  Example: "19:57:22"
  Element: only uses names and attributes of the individual terms
  
%y is replaced by the pattern for the year without century, as a
   decimal number [00,99].
  Pattern: "[0-9][0-9]"
  Example: "89", "01"
  Element name: "year"
  Element attributes: "type" = "short"
  Note: This is the same as %(YY).

%Y is replaced by the pattern for the year, including the century, as a
   decimal number.
  Pattern: "[0-9][0-9][0-9][0-9]"
  Example: "1610", "2002"
  Element name: "year"
  Element attributes: "type" = "long"
  Note: This is the same as %(YYYY).
  
%z is replaced by the pattern for the time-zone as hour offset from GMT.
   (This is used when parsing RFC822-conformant dates, as in
     "%a, %d %b %Y %H:%M:%S %z", except that %z does not include the
    pattern for a missing timezone -- should I fix that?).
  Pattern: "[-+][0-9][0-9][0-9][0-9]"
  Example: "-0500"  (for EST), "+0100" (for CET), "+0530" (somewhere in India)
  Element name: "timezone"
  Element attributes: "type" = "RFC822"

%Z is replaced by a pattern for a timezone name or abbreviation.  (It does
    not allow missing timezone field.)
  Pattern: "(GMT([+-][0-9][0-9][0-9][0-9])?|[A-Z][a-zA-Z]*( [A-Z][a-zA-Z]*)*)"
              (is there anything better?)
  Example: "MST", "GMT", "Pacific Standard Time", "GRNLNDST", "MET DST",
           "New Zealand Standard Time", "NZST", "SAST", "GMT+0200", "IDT"
  Element name: "timezone"
  Element attributes: "type" = "name"
  
%% is replaced by the pattern for "%" (which happens to be "%")
  Pattern: "%"
  Example: "%"
  Element: none

 === Martel specific extensions ===

%(Mon) is the same as "%a".
  Pattern: See the definition for "%a"
  Example: "Wed" "FRI"
  Element name: "weekday"
  Element attributes: "type" = "short"

%(Monday) is the same as "%A".
  Pattern: See the definition for "%A"
  Example: "Thursday" "SUNDAY"
  Element name: "weekday"
  Element attributes: "type" = "long"

%(Jan) is the same as "%b".
  Pattern: See the definition for "%b"
  Example: "Feb"
  Element name: "month"
  Element attributes: "type" = "short"

%(January) is the same as "%B".
  Pattern: See the definition for "%B"
  Example: "August", "MAY"
  Element name: "month"
  Element attributes: "type" = "long"

%(second) is the same as "%S".
  Pattern: See the definition for "%S".
  Example: "03", "25"
  Element name: "second"
  Element attributes: "type" = "numeric"

%(minute) is the same as "%M".
  Pattern: See the definition for "%M"
  Example: "00", "38"
  Element name: "minute"
  Element attributes: "type" = "numeric"

%(12-hour) is replaced by the pattern for a 12 hour clock in any of
   the common formats.  (Numeric values from 1 to 12.)
  Pattern: "(0[1-9]|1[012]?|[2-9]| [1-9])"
  Example: "2", "02", " 2", "10"
  Element name: "hour"
  Element attributes: "type" = "12-hour"

%(24-hour) is replaced by the pattern for a 24 hour clock in any
   of the common formats.  (Numeric values from 0 to 23.)
  Pattern: "([01][0-9]?|2[0123]?|[3-9]| [1-9])"
  Example: "9", "09", " 9", "00", "0", " 0", "23"
  Element name: "hour"
  Element attributes: "type" = "24-hour"

%(hour) is replaced by the pattern for any hour in either a
   12-hour or 24-hour clock.
  Pattern: "([01][0-9]?|2[0123]?|[3-9]| [1-9])"
     (this happens to be the same as %(24-hour)
  Example: "9", "09", " 9", "00", "0", " 0", "23"
  Element name: "hour"
  Element attributes: "type" = "any"

%(day) is replaced by the pattern for the day of the month as a decimal
   in any of the common day format
  Pattern: "(0[1-9]|[12][0-9]?|3[01]?|[4-9]| [1-9])"
  Example: "9", "09", " 9", and "31"
  Element name: "day"
  Element attributes: "type" = "numeric"

%(DD) is the same as "%d", which is the pattern for a day of the month
   as a decimal number [01,31].
  Pattern: See the definition for "%d"
  Example: "09", "31"
  Element name: "day"
  Element attributes: "type" = "numeric"

%(month) is replaced by the pattern for the month as a decimal in any
   of the common month formats.
  Pattern: "(0[1-9]|1[012]?|[2-9]| [1-9])"
  Example: "5", "05", " 5", and "12".
  Element name: "month"
  Element attributes: "type" = "numeric"
  Note: See also "%m" and %(MM).

%(MM) is the same as "%m", which is a two-digit month number [01,12]
  Pattern: See the definition for "%m"
  Example: "05", "01", and "12".
  Element name: "month"
  Element attributes: "type" = "numeric"
  Note: See also %(month).

%(YY)
  Pattern: "[0-9][0-9]"
  Example: "10"
  Element name: "year"
  Element attributes: "type" = "short"

%(YYYY)
  Pattern: "[0-9][0-9][0-9][0-9]"
  Example: "1970"
  Element name: "year"
  Element attributes: "type" = "long"

%(year) is replaced by the pattern accepting 2 digit and 4 digit year formats.
  Pattern: "([0-9]{2}([0-9]{2})?)"
  Example: "2008", "97"
  Element name: "year"
  Element attributes: "type" = "any"
  Note: Need to change this before the year 10,000

"""
import string, Martel, Expression

# Letters are case independent.
def _any_case(s):
    t = ""
    for c in s:
        if c in string.letters:
            t = t + "[%s%s]" % (string.upper(c), string.lower(c))
        else:
            t = t + c
    return t

_time_fields = (
    ("a", _any_case("(Mon|Tue|Wed|Thu|Fri|Sat|Sun)"),
          "weekday", {"type": "short"}),
    ("A", _any_case("(Monday|Tuesday|Wednesday|Thursday|Friday|"
                    "Saturday|Sunday)"),
          "weekday", {"type": "long"}),
    ("b", _any_case("(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)"),
          "month", {"type": "short"}),
    ("B", _any_case("(January|February|March|April|May|June|July|August|"
                   "September|October|November|December)"),
          "month", {"type": "long"}),
    ("C", "\d\d",
          "century", {"type": "numeric"}),
    ("d", "(0[1-9]|[12][0-9]|3[01])", # day of month, "01".."31"
          "day", {"type": "numeric"}),
    ("e", "( [1-9]|[12][0-9]|3[01])", # day of month, " 1".."31"
          "day", {"type": "numeric"}),
    ("g", r"\d{2}",   # ISO 8601 century, 2001 is "01"
          "century", {"type": "ISO8601"}),
    ("G", r"\d{4}",   # ISO 8601 year, 2001
          "year", {"type": "ISO8601"}),
    ("H", "([01][0-9]|2[0-3])", # hours, "00".."23"
          "hour", {"type": "24-hour"}),
    ("I", "(0[0-9]|1[012])",    # hours, "01".."12"
          "hour", {"type": "12-hour"}),

     # day of year, "001".."366"
    ("j", "([12][0-9][0-9]|3([012345][0-9]|6[0-6])|0(0[1-9]|[1-9][0-9]))",
          "year_day", {"type": "1"}),
    
    ("k", "( [0-9]|1[0-9]|2[0123])",  # hour, " 0".."23"
          "hour", {"type": "24-hour"}),
    ("l", "( [0-9]|1[012])",  # hour, " 1", " 2", .. "12"
          "hour", {"type": "12-hour"}),
    ("m", "(0[1-9]|1[012])",  # month, "01", "02", .. "12"
          "month", {"type": "numeric"}),
    ("M", "[0-5][0-9]",       # minute, "00", .. "59"
          "minute", {"type": "numeric"}),
    ("n", r"\n", None, None),
    ("p", "([AaPp][Mm])", # AM
          "ampm", {}),
    ("P", "[aApP][mM]",   # am
          "ampm", {}),
    ("s", r"\d+",    # seconds in unix epoch
          "timestamp", {}),
    ("S", "([0-5][0-9]|6[01])",  # second, [00,61]
          "second", {"type": "numeric"}),
    ("t", r"\t", None, None),
    ("u", "[1-7]",      # weekday, [1,7], 1=Monday
          "weekday", {"type": "Monday1"}),
    ("U", "([0-4][0-9]|5[0-3])",  # week number, [00,53], Sunday first day
          "week_number", {"type": "Sunday_count"}),
    ("V", "(0[1-9]|[1-4][0-9]|5[0-3])",  # week number, [01,53], Monday first day, split
          "week_number", {"type": "type_V"}),  # when is this used?
    ("w", "[0-6]",                # weekday, [0,6], 0=Sunday
          "weekday", {"type": "Sunday0"}),
    ("W", "([0-4][0-9]|5[0-3])",  # week number, [00,53], Monday first day, all
          "week_number", {"type": "Monday_count"}),
    ("y", r"\d{2}",  # 2001 is "01"
          "year", {"type": "short"}),
    ("Y", r"\d{4}",  # 2000
          "year", {"type": "long"}),
    ("z", r"[-+]\d{4}",
          "timezone", {"type": "RFC822"}),
    # "MST", "GMT", "Pacific Standard Time", "GRNLNDST", "MET DST",
    # "New Zealand Standard Time", "NZST", "SAST", "GMT+0200", "IDT"
    ("Z", r"(GMT([+-]\d{4})?|[A-Z][a-zA-Z]*( [A-Z][a-zA-Z]*)*)",
          "timezone", {"type": "name"}),
    ("%", "%", None, None),

    # These are dependent on other fields
    ("D", "%m/%d/%y", None, None),   # 06/29/01
    ("F", "%Y-%m-%d", None, None),   # 2001-06-29
    ("h", "%b", None, None),         # Jan, Feb, ... Dec
    ("r", "%I:%M:%S %p", None, None),# "05:07:50 AM"
    ("R", "%H:%M", None, None),      # "05:07", "19:57"
    ("T", "%H:%M:%S", None, None),   # "01:00:49"
    ("x", "%D", None, None),      # no locale, 09/23/01
    ("X", "%T", None, None),      # no locale, 00:59:18
    ("c", "%a %b %e %T %Y", "date", {}),  # "Wed Dec  2 19:57:22 2001"

    # I made these up
    ("Mon", "%a", None, None),
    ("Monday", "%A", None, None),
    ("Jan", "%b", None, None),
    ("January", "%B", None, None),
    ("second", "%S", None, None),
    ("minute", "%M", None, None),
    ("12-hour", r"(0[1-9]|1[012]?|[2-9]| [1-9])",
           "hour", {"type": "12-hour"}),
    ("24-hour", r"([01][0-9]?|2[0123]?|[3-9]| [0-9])",
           "hour", {"type": "24-hour"}),
    ("hour",   r"([01][0-9]?|2[0123]?|[3-9]| [0-9])",
           "hour", {"type": "any"}),
    ("day", r"(0[1-9]|[12][0-9]?|3[01]?|[4-9]| [1-9])",
              "day", {"type": "numeric"}),
    ("DD", "%d", None, None),
    ("month", r"(0[1-9]|1[012]?|[2-9]| [1-9])", "month", {"type": "numeric"}),
    ("MM", "%m", None, None),
    ("YY", r"[0-9]{2}", "year", {"type": "short"}),
    ("YYYY", r"[0-9]{4}", "year", {"type": "long"}),
    ("year", r"([0-9]{2}([0-9]{2})?)", "year", {"type": "any"}),
)
_time_table = {}
for spec, pat, tag, attrs in _time_fields:
    _time_table[spec] = (pat, tag, attrs)
for v in _time_table.values():
    v = v[0]
    assert (v[0] == '(' and v[-1] == ')') or '|' not in v, v

def make_pattern(format, tag_format = "%s"):
    """format, tag_format = "%s" -> regular expression pattern string

    Turn the given time format string into the corresponding regular
    expression string.  A format term may contain a Group name and attribute
    information.  If present, the group name is %'ed with the
    tag_format to produce the tag name to use.  Use None to specify
    that named groups should not be used.

    >>> from Martel import Time
    >>> print Time.make_pattern("%m-%Y)", "created-%s")
    (?P<created-month?type=numeric>(0[1-9]|1[012]))\\-(?P<created-year?type=long>\\d{4})\\)
    >>>

    See the Time module docstring for more information.
    
    """
    return _parse_time(format, tag_format,
                       text_to_result = Expression.escape,
                       group_to_result = Expression._make_group_pattern,
                       re_to_result = lambda x: x,
                       t = "")

def make_expression(format, tag_format = "%s"):
    """format, tag_format = "%s" -> Martel Expresion

    Turn the given time format string into the corresponding Martel
    Expression.  A format term may contain a Group name and attribute
    information.  If present, the group name is %'ed with the
    tag_format to produce the tag name to use.  Use None to specify
    that named groups should not be used.

    >>> from Martel import Time
    >>> from xml.sax import saxutils
    >>> exp = Time.make_expression("%m-%Y\\n", "created-%s")
    >>> parser = exp.make_parser()
    >>> parser.setContentHandler(saxutils.XMLGenerator())
    >>> parser.parseString("05-1921\n")
    <?xml version="1.0" encoding="iso-8859-1"?>
    <created-month type="numeric">05</created-month>-<created-year type="long">1921</created-year>
    >>> 

    See the Time module docstring for more information.
    
    """
    return _parse_time(format, tag_format,
                       text_to_result = Martel.Str,
                       group_to_result = Martel.Group,
                       re_to_result = Martel.Re,
                       t = Martel.NullOp())

def _use_tag_format(tag_format, name):
    if not tag_format:
        return ""
    return tag_format % name

# text_to_result converts an exact string to the result expression type
#
# group_to_result converts (name, subexpression, attrs) into the
#   result expression type (the subexpression is in the correct type)
#
# re_to_result converts a regular expression pattern string into
#   the result expression type
#
# Partial results are "+"ed to 't'.
#
def _parse_time(s, tag_format, text_to_result, group_to_result,
                re_to_result, t):
    initial_t = t
    n = len(s)
    end = 0
    while end < n:
        prev = end
        # Is there another '%' term?
        start = string.find(s, "%", end)
        if start == -1:
            break
        end = start + 1

        if end == n:
            end = prev
            break # ended with '%' as last character; keep it

        c = s[end]
        # Is this a %() escape?
        if c == '(':
            pos = string.find(s, ")", end)
            if pos != -1:
                # We have a %(special) construct
                c = s[end+1:pos]
            else:
                raise TypeError("Found a '%%(' but no matching ')': %s" % \
                                (repr(s[end-1:]),))
            end = pos

        # Get the expansion and, if needed, do the recursion
        pat, name, attrs = _time_table[c]
        if "%" in pat and pat != "%":
            pat = _parse_time(pat, tag_format,
                              text_to_result, group_to_result,
                              re_to_result, initial_t)
        else:
            pat = re_to_result(pat)
        if name is not None:
            fullname = _use_tag_format(tag_format, name)
            if fullname:
                pat = group_to_result(fullname, pat, attrs)

        if prev + 1 > start:
            t = t + pat
        else:
            t = t + text_to_result(s[prev:start]) + pat
        end = end + 1

    if end < n:
        t = t + text_to_result(s[end:])
    return t

