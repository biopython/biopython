# Check if ReportLab is installed.

try:
    import reportlab as r
    del r
except ImportError:
    raise ImportError("Install ReportLab if you want to use Bio.Graphics. You can find ReportLab at http://www.reportlab.org/downloads.html")
