from Bio.Seq import translate, transcribe, back_transcribe
from wxPython import wx

#TODO - import these from the wx library:
ID_APPLY = 101
ID_CLEAR  = 102
ID_EXIT = 103
ID_CLOSE = 104
ID_ABOUT = 105
ID_CODON = 106
ID_TRANSFORM = 107


class ParamsPanel(wx.wxPanel):
    def __init__(self, parent, log):
        wx.wxPanel.__init__(self, parent, -1)
        codon_table_static = wx.wxStaticText(self, -1, 'Codon Tables',
                                          style=wx.wxALIGN_CENTRE)
        lc = wx.wxLayoutConstraints()
        lc.top.SameAs(self, wx.wxTop, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.wxRight, 5)
        codon_table_static.SetConstraints(lc)

        codon_table_lb = wx.wxListBox(self, ID_CODON)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(codon_table_static, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.PercentOf(self, wx.wxHeight, 30)
        lc.right.SameAs(self, wx.wxRight, 5)
        codon_table_lb.SetConstraints(lc)
        self.codon_table_lb = codon_table_lb

        codon_table_lb.Append('Standard')
        codon_table_lb.Append('Vertebrate Mitochondrial')
        codon_table_lb.Append('Yeast Mitochondrial')
        codon_table_lb.Append('Mold Mitochondrial')
        codon_table_lb.Append('Invertebrate Mitochondrial')
        codon_table_lb.Append('Echinoderm Mitochondrial')
        codon_table_lb.Append('Euplotid Nuclear')
        codon_table_lb.Append('Bacterial')
        codon_table_lb.Append('Alternative Yeast Nuclear')
        codon_table_lb.Append('Ascidian Mitochondrial')
        codon_table_lb.Append('Flatworm Mitochondrial')
        codon_table_lb.Append('Blepharisma Macronuclear')
        codon_table_lb.SetSelection(0)

        transform_static = wx.wxStaticText(self, -1, 'Transformation',
                                        style=wx.wxALIGN_CENTRE)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(codon_table_lb, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.wxRight, 5)
        transform_static.SetConstraints(lc)

        transform_lb = wx.wxListBox(self, ID_TRANSFORM)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(transform_static, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.PercentOf(self, wx.wxHeight, 30)
        lc.right.SameAs(self, wx.wxRight, 5)
        transform_lb.SetConstraints(lc)

        transform_lb.Append('Transcribe')
        transform_lb.Append('Translate')
        transform_lb.Append('Back transcribe')
        transform_lb.SetSelection(1)
        self.transform_lb = transform_lb


class SeqPanel(wx.wxPanel):
    def __init__(self, parent, log):
        self.parent = parent
        wx.wxPanel.__init__(self, parent, -1)
        apply_button = wx.wxButton(self, ID_APPLY, "Apply")
        clear_button = wx.wxButton(self, ID_CLEAR, "Clear")
        close_button = wx.wxButton(self, ID_CLOSE, "Close")
        wx.EVT_BUTTON(self, ID_CLOSE, self.OnClose)
        wx.EVT_BUTTON(self, ID_APPLY, self.OnApply)
        wx.EVT_BUTTON(self, ID_CLEAR, self.OnClear)

        lc = wx.wxLayoutConstraints()
        lc.bottom.SameAs(self, wx.wxBottom, 10)
        lc.left.SameAs(self, wx.wxLeft, 10)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.wxWidth, 25)
        apply_button.SetConstraints(lc)

        lc = wx.wxLayoutConstraints()
        lc.bottom.SameAs(self, wx.wxBottom, 10)
        lc.left.RightOf(apply_button, 5)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.wxWidth, 25)
        clear_button.SetConstraints(lc)

        lc = wx.wxLayoutConstraints()
        lc.bottom.SameAs(self, wx.wxBottom, 10)
        lc.left.RightOf(clear_button, 5)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.wxWidth, 25)
        close_button.SetConstraints(lc)

        src_static = wx.wxStaticText(self, -1, 'Original Sequence',
                                  style=wx.wxALIGN_CENTRE)
        lc = wx.wxLayoutConstraints()
        lc.top.SameAs(self, wx.wxTop, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.wxRight, 5)
        src_static.SetConstraints(lc)

        src_text = wx.wxTextCtrl(self, -1, '', style = wx.wxTE_MULTILINE)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(src_static, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.PercentOf(self, wx.wxHeight, 30)
        lc.right.SameAs(self, wx.wxRight, 5)
        src_text.SetConstraints(lc)
        self.src_text = src_text

        dest_static = wx.wxStaticText(self, -1, 'Transformed Sequence',
                                   style=wx.wxALIGN_CENTRE)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(src_text, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.wxRight, 5)
        dest_static.SetConstraints(lc)

        dest_text = wx.wxTextCtrl(self, -1, '', style = wx.wxTE_MULTILINE)
        lc = wx.wxLayoutConstraints()
        lc.top.Below(dest_static, 5)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.height.PercentOf(self, wx.wxHeight, 30)
        lc.right.SameAs(self, wx.wxRight, 5)
        dest_text.SetConstraints(lc)
        self.dest_text = dest_text


    def OnClose(self, event):
        parent = self.GetParent()
        parent.Destroy()

    def OnApply(self, event):
        codon_table_lb = self.parent.params_panel.codon_table_lb
        selection = codon_table_lb.GetStringSelection()
        print selection
        codon_table = selection[:]
        transform_lb = self.parent.params_panel.transform_lb
        selection = transform_lb.GetStringSelection()
        transform = selection[:]
        print transform
        if(transform == 'Translate'):
            self.translate(codon_table)
        elif(transform == 'Transcribe'):
            self.transcribe()
        elif(transform == 'Back transcribe'):
            self.back_transcribe()

    def OnClear(self, event):
        self.src_text.Clear()
        self.dest_text.Clear()

    def translate(self, codon_table):
        seq = "".join(self.src_text.GetValue().split()) #remove whitespace
        print seq
        self.dest_text.Clear()
        self.dest_text.SetValue(translate(seq, table=codon_table,
                                          to_stop=True))
        
    def transcribe(self):
        seq = "".join(self.src_text.GetValue().split()) #remove whitespace
        print seq
        self.dest_text.Clear()
        self.dest_text.SetValue(transcribe(seq))
                                
    def back_transcribe(self):
        seq = "".join(self.src_text.GetValue().split()) #remove whitespace
        print seq
        self.dest_text.Clear()
        self.dest_text.SetValue(back_transcribe(seq))


class SeqFrame(wx.wxFrame):
    def __init__(self, parent, ID, title):
        wx.wxFrame.__init__(self, parent, ID, title,
                              wx.wxDefaultPosition, wx.wxSize(500, 400))
        self.SetAutoLayout(wx.true)
        self.CreateStatusBar()
        self.SetStatusText("This is the statusbar")
        menu = wx.wxMenu()
        menu.Append(ID_ABOUT, "&About",
                    "More information about this program")
        menu.AppendSeparator()
        menu.Append(ID_EXIT, "E&xit", "Terminate the program")

        menuBar = wx.wxMenuBar()
        menuBar.Append(menu, "&File");
        self.SetMenuBar(menuBar)

        params_panel = ParamsPanel(self, -1)
        lc = wx.wxLayoutConstraints()
        lc.top.SameAs(self, wx.wxTop, 10)
        lc.left.SameAs(self, wx.wxLeft, 5)
        lc.bottom.SameAs(self, wx.wxBottom, 5)
        lc.width.PercentOf(self, wx.wxWidth, 40)
        params_panel.SetConstraints(lc)

        seq_panel = SeqPanel(self, -1)
        lc = wx.wxLayoutConstraints()
        lc.top.SameAs(self, wx.wxTop, 10)
        lc.left.RightOf(params_panel, 5)
        lc.bottom.SameAs(self, wx.wxBottom, 5)
        lc.right.SameAs(self, wx.wxRight)
        seq_panel.SetConstraints(lc)

        self.seq_panel = seq_panel
        self.params_panel = params_panel

        wx.EVT_MENU(self, ID_EXIT, self.exit)

    def exit(self, event):
        self.Close(wx.true)


class MyApp(wx.wxApp):
    def OnInit(self):
        frame = SeqFrame(wx.NULL, -1, "Greetings from Biopython")
        frame.Show(wx.true)
        self.SetTopWindow(frame)
        return wx.true

app = MyApp(0)
app.MainLoop()
