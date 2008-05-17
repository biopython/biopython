# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This code is used to parse XML results in the format specified by
# NCBI's DTD file NCBI_General.mod.dtd (04/10/2008 16:04:22)
# The code is not meant to be used by itself, but is called
# from Bio.Entrez.__init__.py.

error = ""

booleans = (
)

integers = (
    "Date-std_year",		# (%INTEGER;)
    "Date-std_month",		# (%INTEGER;)
    "Date-std_day",		# (%INTEGER;)
    "Date-std_hour",		# (%INTEGER;)
    "Date-std_minute",		# (%INTEGER;)
    "Date-std_second",		# (%INTEGER;)
    "User-field_data_ints_E",	# (%INTEGER;)
    "User-field_data_int",	# (%INTEGER;)
    "User-field_num",		# (%INTEGER;)
    "Int-fuzz_alt_E",		# (%INTEGER;)
    "Int-fuzz_range_max",	# (%INTEGER;)
    "Int-fuzz_p-m",		# (%INTEGER;)
    "Object-id_id",		# (%INTEGER;)
    "Int-fuzz_range_min",	# (%INTEGER;)
    "Int-fuzz_pct",		# (%INTEGER;)
)

strings = (
    "Date_str",			# (#PCDATA)
    "User-field_data_strs_E",	# (#PCDATA)
    "User-field_data_str",	# (#PCDATA)
    "User-object_class",	# (#PCDATA)
    "Name-std_first",		# (#PCDATA)
    "Name-std_middle",		# (#PCDATA)
    "Name-std_full",		# (#PCDATA)
    "Name-std_last",		# (#PCDATA)
    "Person-id_ml",		# (#PCDATA)
    "Object-id_str",		# (#PCDATA)
    "Dbtag_db",			# (#PCDATA)
    "Date-std_season",		# (#PCDATA)
    "Person-id_str",		# (#PCDATA)
    "Person-id_consortium",	# (#PCDATA)
    "Name-std_initials",	# (#PCDATA)
    "Name-std_suffix",		# (#PCDATA)
    "Name-std_title",		# (#PCDATA)
    "Int-fuzz_lim",		# %ENUM;
				# ATTLIST value (unk | gt | lt | tr | tl
				#               | circle | other) #REQUIRED
    "User-field_data_real",	# (%REAL;)
    "User-field_data_bool",	# EMPTY
				# ATTLIST value ( true | false ) #REQUIRED
    "User-field_data_os",	# (%OCTETS;)
    "User-field_data_reals_E",	# (%REAL;)
    "User-field_data_oss_E",	# (%OCTETS;)
)

lists = (
    "User-field_data_reals",	# (User-field_data_reals_E*)
    "User-field_data_fields",	# (User-field*)
    "User-field_data_objects",	# (User-object*)
    "User-field_data_oss",	# (User-field_data_oss_E*)
    "User-field_data_ints",	# (User-field_data_ints_E*)
    "User-field_data_strs",	# (User-field_data_strs_E*)
    "User-object_data",		# (User-field*)
    "Int-fuzz_alt",		# (Int-fuzz_alt_E*)
)

dictionaries = (
    "Date",		# ( Date_str | Date_std)
    "Date_std",		# (Date-std>
    "Date-std",		# (Date-std_year, 
			#  Date-std_month?, 
			#  Date-std_day?, 
			#  Date-std_season?, 
			#  Date-std_hour?, 
			#  Date-std_minute?, 
			#  Date-std_second?)
    "Dbtag",		# ( Dbtag_db, Dbtag_tag)
    "Dbtag_tag",	# (Object-id)
    "Object-id",	# ( Object-id_id | Object-id_str)
    "Person-id",	# ( Person-id_dbtag | 
			#   Person-id_name | 
			#   Person-id_ml | 
			#   Person-id_str | 
			#   Person-id_consortium)
    "Person-id_dbtag",	# (Dbtag)
    "Person-id_name",	# (Name-std)
    "Name-std",		# (Name-std_last, 
			#  Name-std_first?, 
			#  Name-std_middle?, 
			#  Name-std_full?, 
			#  Name-std_initials?, 
			#  Name-std_suffix?, 
			#  Name-std_title?)
    "Int-fuzz",		# (Int-fuzz_p-m | 
			#  Int-fuzz_range | 
			#  Int-fuzz_pct | 
			#  Int-fuzz_lim | 
			#  Int-fuzz_alt)
    "Int-fuzz_range",	# (Int-fuzz_range_max, Int-fuzz_range_min)
    "User-object",	# (User-object_class?, 
			#  User-object_type, 
			#  User-object_data)
    "User-object_type",	# (Object-id)
    "User-field",	# (User-field_label, User-field_num?, User-field_data)
    "User-field_label",	# (Object-id)
    "User-field_data",	# (User-field_data_str | 
			#  User-field_data_int | 
			#  User-field_data_real | 
			#  User-field_data_bool | 
			#  User-field_data_os | 
			#  User-field_data_object | 
			#  User-field_data_strs | 
			#  User-field_data_ints | 
			#  User-field_data_reals | 
			#  User-field_data_oss | 
			#  User-field_data_fields | 
			#  User-field_data_objects)
    "User-field_data_object",	# (User-object)
)

structures = (
)

items = (
)


def startElement(self, name, attrs):
    return

def endElement(self, name):
    self.path = self.path[:-1]
