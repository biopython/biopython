from Martel import Dispatch
from Bio import Search, StdHandler

class BuildSearch(Dispatch.Dispatcher):
    def __init__(self):
        Dispatch.Dispatcher.__init__(self,
                                     prefix = "bioformat:",
                                     remap = {
            "bioformat:hit": "bioformat:hit_description_block",
            }
                                     )
        
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
        self.document = None
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

    def add_hsp(self, hsp_values, hsp_handler, strands, frames):
        self.hsps.append(Search.HSP(hsp_handler.query_seq,
                                    hsp_handler.homology_seq,
                                    hsp_handler.subject_seq,

                                    # XXX strand and frame!
                                    (hsp_handler.query_start_loc,
                                     hsp_handler.query_end_loc),
                                    (hsp_handler.subject_start_loc,
                                     hsp_handler.subject_end_loc),

                                    hsp_handler.query_name,
                                    hsp_handler.subject_name,

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

StdHandler.add_text_block_handler(BuildSearch, "hit_description",
                                  "join-description", "join|fixspaces",
                                  "hit_description")

StdHandler.add_int_handler(BuildSearch, "hit_length", "hit_length")
                                  
make_builder = BuildSearch
