from wxPython.wx import *

ID_ABOUT = 101
ID_EXIT  = 102

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

        params_panel = wxPanel(self, -1)
	lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 10 )
        lc.left.SameAs( self, wxLeft, 5 )
        lc.bottom.SameAs( self, wxBottom, 5 )
        lc.width.PercentOf( self, wxWidth, 40 )
        params_panel.SetConstraints( lc )

        seq_panel = wxPanel(self, -1)
	lc = wxLayoutConstraints()
        lc.top.SameAs( self, wxTop, 10 )
        lc.left.RightOf( params_panel, 5 )
        lc.bottom.SameAs( self, wxBottom, 5 )
        lc.right.SameAs( self, wxRight )
        seq_panel.SetConstraints( lc )



        translate_button = wxButton( self, -1, "Apply" )
        transcribe_button = wxButton( self, -1, "Clear" )
        reverse_button = wxButton( self, -1, "Close" )

	lc = wxLayoutConstraints()
        lc.bottom.SameAs( seq_panel, wxBottom, 10 )
        lc.left.SameAs( seq_panel, wxLeft, 10 )
        lc.height.AsIs( )
        lc.width.PercentOf( seq_panel, wxWidth, 25 )
        translate_button.SetConstraints( lc )

        lc = wxLayoutConstraints()
        lc.bottom.SameAs( seq_panel, wxBottom, 10 )
        lc.left.RightOf( translate_button, 5 )
        lc.height.AsIs()
        lc.width.PercentOf( seq_panel, wxWidth, 25 )
        transcribe_button.SetConstraints( lc )

        lc = wxLayoutConstraints()
        lc.bottom.SameAs( seq_panel, wxBottom, 10 )
        lc.left.RightOf( transcribe_button, 5 )
        lc.height.AsIs()
        lc.width.PercentOf( seq_panel, wxWidth, 25 )
        reverse_button.SetConstraints( lc )

        src_static = wxStaticText( seq_panel, -1, 'Original Sequence', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.SameAs( seq_panel, wxTop, 5 )
        lc.left.SameAs( seq_panel, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( seq_panel, wxRight, 5 )
        src_static.SetConstraints( lc )

        src_text = wxTextCtrl( seq_panel, -1, '', style = wxTE_MULTILINE )
        lc = wxLayoutConstraints()
        lc.top.Below( src_static, 5 )
        lc.left.SameAs( seq_panel, wxLeft, 5 )
        lc.height.PercentOf( seq_panel, wxHeight, 30 )
        lc.right.SameAs( seq_panel, wxRight, 5 )
        src_text.SetConstraints( lc )

        dest_static = wxStaticText( seq_panel, -1, 'Transformed Sequence', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.Below( src_text, 5 )
        lc.left.SameAs( seq_panel, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( seq_panel, wxRight, 5 )
        dest_static.SetConstraints( lc )

        dest_text = wxTextCtrl( seq_panel, -1, '', style = wxTE_MULTILINE )
        lc = wxLayoutConstraints()
        lc.top.Below( dest_static, 5 )
        lc.left.SameAs( seq_panel, wxLeft, 5 )
        lc.height.PercentOf( seq_panel, wxHeight, 30 )
        lc.right.SameAs( seq_panel, wxRight, 5 )
        dest_text.SetConstraints( lc )

        codon_table_static = wxStaticText( params_panel, -1, 'Codon Tables', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.SameAs( params_panel, wxTop, 5 )
        lc.left.SameAs( params_panel, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( params_panel, wxRight, 5 )
        codon_table_static.SetConstraints( lc )

        codon_table_lb = wxListBox( params_panel, -1 )
        lc = wxLayoutConstraints()
        lc.top.Below( codon_table_static, 5 )
        lc.left.SameAs( params_panel, wxLeft, 5 )
        lc.height.PercentOf( params_panel, wxHeight, 30 )
        lc.right.SameAs( params_panel, wxRight, 5 )
        codon_table_lb.SetConstraints( lc )

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

        transform_static = wxStaticText( params_panel, -1, 'Transformation', \
	    style = wxALIGN_CENTRE )
        lc = wxLayoutConstraints()
        lc.top.Below( codon_table_lb, 5 )
        lc.left.SameAs( params_panel, wxLeft, 5 )
        lc.height.AsIs()
        lc.right.SameAs( params_panel, wxRight, 5 )
        transform_static.SetConstraints( lc )

        transform_lb = wxListBox( params_panel, -1 )
        lc = wxLayoutConstraints()
        lc.top.Below( transform_static, 5 )
        lc.left.SameAs( params_panel, wxLeft, 5 )
        lc.height.PercentOf( params_panel, wxHeight, 30 )
        lc.right.SameAs( params_panel, wxRight, 5 )
        transform_lb.SetConstraints( lc )

        transform_lb.Append( 'Transcribe' )
        transform_lb.Append( 'Translate' )
        transform_lb.Append( 'Back translate' )




class MyApp(wxApp):
    def OnInit(self):
        frame = SeqFrame(NULL, -1, "Greetings from biopython")
        frame.Show(true)
        self.SetTopWindow(frame)
        return true

app = MyApp(0)
app.MainLoop()
