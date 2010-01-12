from Bio.Seq import translate, transcribe, back_transcribe
import wx

ID_APPLY = 101
ID_CLEAR  = 102
ID_EXIT = 103
ID_CLOSE = 104
ID_ABOUT = 105
ID_CODON = 106
ID_TRANSFORM = 107


class ParamsPanel(wx.Panel):
    def __init__(self, parent, log):
        wx.Panel.__init__(self, parent, -1)
        codon_table_static = wx.StaticText(self, -1, 'Codon Tables',
                                           style=wx.ALIGN_CENTER)
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.Right, 5)
        codon_table_static.SetConstraints(lc)

        codon_table_lb = wx.ListBox(self, ID_CODON)
        lc = wx.LayoutConstraints()
        lc.top.Below(codon_table_static, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.PercentOf(self, wx.Height, 30)
        lc.right.SameAs(self, wx.Right, 5)
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

        transform_static = wx.StaticText(self, -1, 'Transformation',
                                         style=wx.ALIGN_CENTER)
        lc = wx.LayoutConstraints()
        lc.top.Below(codon_table_lb, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.Right, 5)
        transform_static.SetConstraints(lc)

        transform_lb = wx.ListBox(self, ID_TRANSFORM)
        lc = wx.LayoutConstraints()
        lc.top.Below(transform_static, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.PercentOf(self, wx.Height, 30)
        lc.right.SameAs(self, wx.Right, 5)
        transform_lb.SetConstraints(lc)

        transform_lb.Append('Transcribe')
        transform_lb.Append('Translate')
        transform_lb.Append('Back transcribe')
        transform_lb.SetSelection(1)
        self.transform_lb = transform_lb


class SeqPanel(wx.Panel):
    def __init__(self, parent, log):
        self.parent = parent
        wx.Panel.__init__(self, parent, -1)
        apply_button = wx.Button(self, ID_APPLY, "Apply")
        clear_button = wx.Button(self, ID_CLEAR, "Clear")
        close_button = wx.Button(self, ID_CLOSE, "Close")
        wx.EVT_BUTTON(self, ID_CLOSE, self.OnClose)
        wx.EVT_BUTTON(self, ID_APPLY, self.OnApply)
        wx.EVT_BUTTON(self, ID_CLEAR, self.OnClear)

        lc = wx.LayoutConstraints()
        lc.bottom.SameAs(self, wx.Bottom, 10)
        lc.left.SameAs(self, wx.Left, 10)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.Width, 25)
        apply_button.SetConstraints(lc)

        lc = wx.LayoutConstraints()
        lc.bottom.SameAs(self, wx.Bottom, 10)
        lc.left.RightOf(apply_button, 5)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.Width, 25)
        clear_button.SetConstraints(lc)

        lc = wx.LayoutConstraints()
        lc.bottom.SameAs(self, wx.Bottom, 10)
        lc.left.RightOf(clear_button, 5)
        lc.height.AsIs()
        lc.width.PercentOf(self, wx.Width, 25)
        close_button.SetConstraints(lc)

        src_static = wx.StaticText(self, -1, 'Original Sequence',
                                   style=wx.ALIGN_CENTER)
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.Right, 5)
        src_static.SetConstraints(lc)

        src_text = wx.TextCtrl(self, -1, '', style=wx.TE_MULTILINE)
        lc = wx.LayoutConstraints()
        lc.top.Below(src_static, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.PercentOf(self, wx.Height, 30)
        lc.right.SameAs(self, wx.Right, 5)
        src_text.SetConstraints(lc)
        self.src_text = src_text

        dest_static = wx.StaticText(self, -1, 'Transformed Sequence',
                                    style=wx.ALIGN_CENTER)
        lc = wx.LayoutConstraints()
        lc.top.Below(src_text, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.AsIs()
        lc.right.SameAs(self, wx.Right, 5)
        dest_static.SetConstraints(lc)

        dest_text = wx.TextCtrl(self, -1, '', style=wx.TE_MULTILINE)
        lc = wx.LayoutConstraints()
        lc.top.Below(dest_static, 5)
        lc.left.SameAs(self, wx.Left, 5)
        lc.height.PercentOf(self, wx.Height, 30)
        lc.right.SameAs(self, wx.Right, 5)
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


class SeqFrame(wx.Frame):
    def __init__(self, parent, ID, title):
        wx.Frame.__init__(self, parent, ID, title,
                          wx.DefaultPosition, wx.Size(500, 400))
        self.SetAutoLayout(True)
        self.CreateStatusBar()
        self.SetStatusText("This is the statusbar")
        menu = wx.Menu()
        menu.Append(ID_ABOUT, "&About",
                    "More information about this program")
        menu.AppendSeparator()
        menu.Append(ID_EXIT, "E&xit", "Terminate the program")

        menuBar = wx.MenuBar()
        menuBar.Append(menu, "&File");
        self.SetMenuBar(menuBar)

        params_panel = ParamsPanel(self, -1)
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 10)
        lc.left.SameAs(self, wx.Left, 5)
        lc.bottom.SameAs(self, wx.Bottom, 5)
        lc.width.PercentOf(self, wx.Width, 40)
        params_panel.SetConstraints(lc)

        seq_panel = SeqPanel(self, -1)
        lc = wx.LayoutConstraints()
        lc.top.SameAs(self, wx.Top, 10)
        lc.left.RightOf(params_panel, 5)
        lc.bottom.SameAs(self, wx.Bottom, 5)
        lc.right.SameAs(self, wx.Right)
        seq_panel.SetConstraints(lc)

        self.seq_panel = seq_panel
        self.params_panel = params_panel

        wx.EVT_MENU(self, ID_EXIT, self.exit)

    def exit(self, event):
        self.Close(True)


class MyApp(wx.App):
    def OnInit(self):
        frame = SeqFrame(None, -1, "Greetings from Biopython")
        frame.Show(True)
        self.SetTopWindow(frame)
        return True

app = MyApp(0)
app.MainLoop()
