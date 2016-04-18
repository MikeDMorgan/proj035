import CGATPipelines.Pipeline as P
import CGAT.Experiment as E
import rpy2.robjects as ro
from rpy2.robjects import r as R
import pandas as pd
import pandas.rpy.common as com
import os
import re
import PipelineProject035 as P35

print "hello world"

new_df = pd.read_table("design.tsv",
                       sep="\t", index_col=0, header=0)

print new_df.head()
