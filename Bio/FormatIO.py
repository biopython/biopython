#Would like to have just issued a deprecation warning, and removed this
#module later.  However, due to the FormatIO code in Bio/SeqRecord.py the
#deprecation warning would be triggered whenever someone used the SeqRecord.
raise ImportError, "Bio.FormatIO has been removed.  Please try Bio.SeqIO instead"
