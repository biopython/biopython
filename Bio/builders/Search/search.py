from Bio import Search, Dispatch, StdHandler

class BuildSearch(Dispatch.Dispatcher):
    def __init__(self):
        Dispatch.Dispatcher.__init__(self,
                                     prefix = "bioformat:",
                                     remap = {
            "bioformat:hit": "bioformat:hit_description_block",
            })
        
        self.acquire(StdHandler.TextBlock(self, "hit_description"),
                     "hit_description")
        self.acquire(StdHandler.Int(self, "hit_length"), "hit_length")
                     
        
        self.acquire(StdHandler.Handle_hsp(self.add_hsp))
        self.acquire(StdHandler.Handle_search_header(self.add_header))
        self.acquire(StdHandler.Handle_search_table(self.add_table))
        self.acquire(StdHandler.Handle_search_info(self.add_stats))

    def start_(self, name, attrs):
        self.algorithm = None
        self.query = None
        self.database = None
        self.table = []
        self.hits = []
        self.parameters = {}
        self.statistics = {}
        
    def end_(self, name):
        self.document = Search.Search(
            self.algorithm,
            self.query,
            self.database,
            self.table,
            self.hits,
            self.parameters,
            self.statistics)

    def start_hit(self, name, attrs):
        self.hit_description = None
        self.hit_length = None
        self.hsps = []

    def end_hit(self, name):
        self.hits.append(Search.Hit("XXX SPAM", self.hit_description,
                                    "XXX EGGS", self.hit_length,
                                    self.algorithm, self.hsps))

    def add_hsp(self, hsp_values, names, seqs, starts, ends, strands, frames):
        
        self.hsps.append(Search.HSP(seqs["query"],
                                    seqs["homology"],
                                    seqs["subject"],

                                    # XXX strand and frame!
                                    (starts["query"], ends["query"]),
                                    (starts["subject"], ends["subject"]),

                                    names["query"],
                                    names["subject"],

                                    self.algorithm,

                                    hsp_values))

    def add_table(self, table):
        self.table = [Search.TableInfo(*x) for x in table]

    def add_header(self, info):
        self.algorithm = Search.Algorithm(info["appname"],
                                          info["appversion"])
        self.database = Search.Database(info["dbname"],
                                        info["db_num_letters"],
                                        info["db_num_sequences"])
        self.query = Search.Query("XXX spam", "XXX eggs",
                                  info["query_description"],
                                  info["query_size"])

    def add_stats(self, parameters, statistics):
        self.parameters = parameters
        self.statistics = statistics

make_builder = BuildSearch
