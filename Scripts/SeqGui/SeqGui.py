import string
from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from wxPython.wx import *
from Bio import Translate
from Bio import Transcribe

ID_APPLY = 101
ID_CLEAR  = 102
ID_EXIT = 103
ID_CLOSE = 104
ID_ABOUT = 105
ID_CODON = 106
ID_TRANSFORM = 107


class ParamsPanel( wxPanel ):
    def __init__(self, parent, log):
        wxPanel.__init__(self, parent, -1)
        codon_table_static = wxStaticText( self, -1, 'Codon Tables', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( self, wxRight, 5 )
        codon_table_static.SetConstraints( lc )

        codon_table_lb = wxListBox( self, ID_CODON )
        lc = wxLayoutConstraints()
        lc.top.Below( codon_table_static, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.PercentOf( self, wxHeight, 30 )
        lc.right.SameAs( self, wxRight, 5 )
        codon_table_lb.SetConstraints( lc )
        self.codon_table_lb = codon_table_lb

        codon_table_lb.Append( 'Standard' )
        codon_table_lb.Append( 'Vertebrate Mitochondrial' )
        codon_table_lb.Append( 'Yeast Mitochondrial' )
        codon_table_lb.Append( 'Mold Mitochondrial' )
        codon_table_lb.Append( 'Invertebrate Mitochondrial' )
        codon_table_lb.Append( 'Echinoderm Mitochondrial' )
        codon_table_lb.Append( 'Euplotid Nuclear' )
        codon_table_lb.Append( 'Bacterial' )
        codon_table_lb.Append( 'Alternative Yeast Nuclear' )
        codon_table_lb.Append( 'Ascidian Mitochondrial' )
        codon_table_lb.Append( 'Flatworm Mitochondrial' )
        codon_table_lb.Append( 'Blepharisma Macronuclear' )
        codon_table_lb.SetSelection( 0 )

        transform_static = wxStaticText( self, -1, 'Transformation', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.Below( codon_table_lb, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( self, wxRight, 5 )
        transform_static.SetConstraints( lc )

        transform_lb = wxListBox( self, ID_TRANSFORM )
        lc = wxLayoutConstraints()
        lc.top.Below( transform_static, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.PercentOf( self, wxHeight, 30 )
        lc.right.SameAs( self, wxRight, 5 )
        transform_lb.SetConstraints( lc )

        transform_lb.Append( 'Transcribe' )
        transform_lb.Append( 'Translate' )
        transform_lb.Append( 'Back translate' )
        transform_lb.Append( 'Back transcribe' )
        transform_lb.SetSelection( 1 )
        self.transform_lb = transform_lb


class SeqPanel( wxPanel ):
    def __init__(self, parent, log):
        self.parent = parent
        wxPanel.__init__(self, parent, -1)
        apply_button = wxButton( self, ID_APPLY, "Apply" )
        clear_button = wxButton( self, ID_CLEAR, "Clear" )
        close_button = wxButton( self, ID_CLOSE, "Close" )
        EVT_BUTTON( self, ID_CLOSE, self.OnClose )
        EVT_BUTTON( self, ID_APPLY, self.OnApply )
        EVT_BUTTON( self, ID_CLEAR, self.OnClear )

	lc = wxLayoutConstraints()
        lc.bottom.SameAs( self, wxBottom, 10 )
        lc.left.SameAs( self, wxLeft, 10 )
        lc.height.AsIs( )
        lc.width.PercentOf( self, wxWidth, 25 )
        apply_button.SetConstraints( lc )

        lc = wxLayoutConstraints()
        lc.bottom.SameAs( self, wxBottom, 10 )
        lc.left.RightOf( apply_button, 5 )
        lc.height.AsIs()
        lc.width.PercentOf( self, wxWidth, 25 )
        clear_button.SetConstraints( lc )

        lc = wxLayoutConstraints()
        lc.bottom.SameAs( self, wxBottom, 10 )
        lc.left.RightOf( clear_button, 5 )
        lc.height.AsIs()
        lc.width.PercentOf( self, wxWidth, 25 )
        close_button.SetConstraints( lc )

        src_static = wxStaticText( self, -1, 'Original Sequence', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( self, wxRight, 5 )
        src_static.SetConstraints( lc )

        src_text = wxTextCtrl( self, -1, '', style = wxTE_MULTILINE )
        lc = wxLayoutConstraints()
        lc.top.Below( src_static, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.PercentOf( self, wxHeight, 30 )
        lc.right.SameAs( self, wxRight, 5 )
        src_text.SetConstraints( lc )
        self.src_text = src_text

        dest_static = wxStaticText( self, -1, 'Transformed Sequence', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.Below( src_text, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( self, wxRight, 5 )
        dest_static.SetConstraints( lc )

        dest_text = wxTextCtrl( self, -1, '', style = wxTE_MULTILINE )
        lc = wxLayoutConstraints()
        lc.top.Below( dest_static, 5 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.height.PercentOf( self, wxHeight, 30 )
        lc.right.SameAs( self, wxRight, 5 )
        dest_text.SetConstraints( lc )
        self.dest_text = dest_text


    def OnClose( self, event ):
        parent = self.GetParent()
        parent.Destroy()

    def OnApply( self, event ):
        codon_table_lb = self.parent.params_panel.codon_table_lb
        selection = codon_table_lb.GetStringSelection()
        print selection
        codon_table = selection[:]
        transform_lb = self.parent.params_panel.transform_lb
        selection = transform_lb.GetStringSelection()
        transform = selection[:]
        print transform
        if( transform == 'Translate' ):
            self.translate( codon_table )
        elif( transform == 'Back translate' ):
            self.back_translate( codon_table )
        elif( transform == 'Transcribe' ):
            self.transcribe()
        elif( transform == 'Back transcribe' ):
            self.back_transcribe()

    def OnClear( self, event ):
        self.src_text.Clear()
        self.dest_text.Clear()

    def translate( self, codon_table ):
        trans = Translate.unambiguous_dna_by_name[ codon_table ]
        text = self.src_text.GetValue()
        seq = text[:]
        seq = string.join( string.split( seq ) )
        dna = Seq.Seq( seq, IUPAC.unambiguous_dna )
        print dna
        protein = trans.translate_to_stop( dna )
        self.dest_text.Clear()
        self.dest_text.SetValue( protein.tostring() )

    def back_translate( self, codon_table ):
        trans = Translate.unambiguous_dna_by_name[ codon_table ]
        text = self.src_text.GetValue()
        seq = text[:]
        seq = string.join( string.split( seq ) )
        protein = Seq.Seq( seq, IUPAC.unambiguous_dna )
        print protein
        dna = trans.back_translate( protein )
        self.dest_text.Clear()
        self.dest_text.SetValue( dna.tostring() )

    def transcribe( self ):
        trans = Transcribe.unambiguous_transcriber
        text = self.src_text.GetValue()
        seq = text[:]
        seq = string.join( string.split( seq ) )
        dna = Seq.Seq( seq, IUPAC.unambiguous_dna )
        print dna
        rna = trans.transcribe( dna )
        self.dest_text.Clear()
        self.dest_text.SetValue( rna.tostring() )

    def back_transcribe( self ):
        trans = Transcribe.unambiguous_transcriber
        text = self.src_text.GetValue()
        seq = text[:]
        seq = string.join( string.split( seq ) )
        rna = Seq.Seq( seq, IUPAC.unambiguous_rna )
        print rna
        dna = trans.back_transcribe( rna )
        self.dest_text.Clear()
        self.dest_text.SetValue( dna.tostring() )





class SeqFrame(wxFrame):
    def __init__(self, parent, ID, title):
        wxFrame.__init__(self, parent, ID, title,
                              wxDefaultPosition, wxSize(500, 400))
        self.SetAutoLayout( true )
	self.CreateStatusBar()
        self.SetStatusText("This is the statusbar")
        menu = wxMenu()
        menu.Append(ID_ABOUT, "&About",
                         "More information about this program")
        menu.AppendSeparator()
        menu.Append(ID_EXIT, "E&xit", "Terminate the program")

        menuBar = wxMenuBar()
        menuBar.Append(menu, "&File");
        self.SetMenuBar(menuBar)

        params_panel = ParamsPanel(self, -1)
	lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 10 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.bottom.SameAs( self, wxBottom, 5 )
        lc.width.PercentOf( self, wxWidth, 40 )
        params_panel.SetConstraints( lc )

        seq_panel = SeqPanel(self, -1)
	lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 10 )
        lc.left.RightOf( params_panel, 5 )
        lc.bottom.SameAs( self, wxBottom, 5 )
        lc.right.SameAs( self, wxRight )
        seq_panel.SetConstraints( lc )

        self.seq_panel = seq_panel
        self.params_panel = params_panel



        EVT_MENU( self, ID_EXIT, self.exit )

    def exit( self, event ):
        self.Close( true )




class MyApp(wxApp):
    def OnInit(self):
        frame = SeqFrame(NULL, -1, "Greetings from biopython")
        frame.Show(true)
        self.SetTopWindow(frame)
        return true

app = MyApp(0)
app.MainLoop()
