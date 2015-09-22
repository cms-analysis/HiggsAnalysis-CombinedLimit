class PhysOpts(object):
  def __init__(self):
    self.opts = set()
  def add(self, name, default=0.):
    setattr(self, name, default)
    self.opts.add(name)
  def parse(self, opt):
    key,value =tuple(opt.split('='))
    if key in self.opts:
      t = type(getattr(self, key))
      if isinstance(getattr(self, key), bool):
         vals = {'True' : True, 'False' : False}
         if value in vals:
            setattr(self, key, vals[value])
         else:
            raise ValueError(
               'Model option %s was declared '
               'as %s but got %s!' % (key, t, value)
               )
      else:
         try:
            setattr(self, key, t(value))
         except ValueError:
            raise ValueError(
               'Model option %s was declared '
               'as %s but got %s!' % (key, t, value)
               )
    else:
      raise ValueError(
         "Model option %s not recognised!"
         " Available options are: %s" %(key, ','.join(self.opts))
         )
