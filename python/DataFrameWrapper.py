from __future__ import absolute_import

import numpy as np
import pandas as pd
from six.moves import range

import ROOT


class DataFrameWrapper(object):
    th1_class = "TH1F"

    def __init__(self, path, ext, load=True):
        """
        path: file system path to the dataframe on disk
        ext: file extension used to interpret the file based on pandas IO
        load: skip loading the dataframe if False
        """
        self.path = path
        self.ext = ext

        self.read_args = []
        self.read_kwargs = {}
        if load:
            self.df = self.load_dataframe()

    def load_dataframe(self):
        """
        Use pandas IO tools to load a dataframe from any of the following:
        ["csv", "json", "html", "pkl", "xlsx", "h5", "parquet"]
        """
        # Index all columns apart from the last 2 (taken as sum_w and sum_ww)
        # for csv, json, html and xlsx
        if self.ext == ".csv":
            df = pd.read_csv(self.path, *self.read_args, **self.read_kwargs)
            return df.set_index(df.columns.tolist()[:-2])
        elif self.ext == ".json":
            df = pd.read_json(self.path, *self.read_args, **self.read_kwargs)
            return df.set_index(df.columns.tolist()[:-2])
        elif self.ext == ".html":
            df = pd.read_html(self.path, *self.read_args, **self.read_kwargs)
            return df.set_index(df.columns.tolist()[:-2])
        elif self.ext == ".pkl":
            return pd.read_pickle(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".xlsx":
            if ":" in self.path:
                filepath, sheetname = self.path.split(":")
            else:
                filepath, sheetname = self.path, 0

            # read in columns first
            cols = pd.read_excel(self.path, sheetname, *self.read_args, **self.read_kwargs).columns.tolist()
            return pd.read_excel(self.path, sheetname, index_col=list(range(len(cols) - 2)), *self.read_args, **self.read_kwargs)
        elif self.ext == ".h5":
            filepath, internalpath = self.path.split(":")
            return pd.read_hdf(filepath, internalpath, *self.read_args, **self.read_kwargs)
        elif self.ext == ".parquet":
            return pd.read_parquet(self.path, *self.read_args, **self.read_kwargs)
        return None

    def Get(self, object_name):
        """
        Mimic ROOT file Get function to return a ROOT.TH1.

        object_name: dataframe index/columns selection split by ":" (used in .loc)
            e.g. "125:bin1:signal:sigmaUp,event_count:event_variance" to select
            the index (125, bin1, signal, sigmaUp) and columns (event_count,
            event_variance)
        """
        if not hasattr(self, "df"):
            raise AttributeError("Dataframe has not been loaded")

        # index_labels, column_labels = object_name.split(",")
        # Not necessary to include ":nominal,sum_w:sum_ww in every line
        # check if its not there and if so, add these default names of columns
        try:
            index_labels, column_labels = object_name.split(",")
        except ValueError:
            index_labels = object_name
            column_labels = ""
        if len(column_labels) == 0:
            column_labels = "sum_w:sum_ww"
        if len(index_labels.split(":")) < 3:
            index_labels += ":nominal"

        # Try to cast index_labels into self.df.index dtypes. Users can only
        # input index_labels as a string, but the dataframe might have other
        # dtypes (e.g. int, float, ...)
        selection = [
            np.array([index_labels.split(":")[idx]]).astype(self.df.index.get_level_values(idx).dtype)[0] for idx in range(len(index_labels.split(":")))
        ]

        df_hist = self.df.loc[tuple(selection), column_labels.split(":")]

        # index name used for th1 name
        df_hist.index.names = [index_labels.replace(":", "_")]
        return self.convert_to_th1(df_hist, self.th1_class)

    @staticmethod
    def convert_to_th1(df, th1_class):
        """
        Receive a dataframe and convert it to a TH1. Index is taken as the
        binning for labelling. Overflowbin should be ignored!
        """
        name = df.index.names[0]
        nbins = df.shape[0]
        th1 = getattr(ROOT, th1_class)(name, name, nbins, 0.0, float(nbins))
        for i in range(nbins):
            sum_w, sum_ww = df.iloc[i]
            th1.SetBinContent(i + 1, sum_w)
            th1.SetBinError(i + 1, np.sqrt(sum_ww))
        return th1
