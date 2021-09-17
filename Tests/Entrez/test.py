from Bio import Entrez

Entrez.email = "mjldehoon@yahoo.com"

# stream = Entrez.efetch(db='pmc', id="8435807")
# data = stream.read()
# stream.close()
# handle =  open("efetch_pmc.xml", 'wb')
# handle.write(data)
# handle.close()

handle = open("efetch_pmc.xml", "rb")
record = Entrez.read(handle)
handle.close()
