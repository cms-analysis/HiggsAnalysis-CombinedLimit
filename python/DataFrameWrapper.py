import pandas as pd

class DataFrameWrapper(object):
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
        if self.ext == ".csv":
            return pd.read_csv(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".json":
            return pd.read_json(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".html":
            return pd.read_html(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".pkl":
            return pd.read_pickle(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".xlsx":
            return pd.read_excel(self.path, *self.read_args, **self.read_kwargs)
        elif self.ext == ".h5":
            filepath, internalpath = self.path.split(":")
            return pd.read_hdf(filepath, internalpath, *self.read_args, **self.read_kwargs)
        elif self.ext == ".parquet":
            return pd.read_parquet(self.path, *self.read_args, **self.read_kwargs)
        return None

    def Get(self, object_name):
        """
        Mimic ROOT file Get function to return a ROOT.TH1

        object_name: dataframe index selection split by ":" (used in .loc)
        """
        raise NotImplementedError()
