# -*- coding: utf-8 -*-
"""
LS-DYNA keyfile parser
Tomasz Koziara 2014-08-22
"""
import sys, os, re

def Keyfile (path):
  """ parse LS-DYNA keyfile """
  
  # supported keywords
  keywords = set (['*KEYWORD', '*DEFINE_CURVE', '*PART', '*PART_INERTIA', '*MAT_RIGID', '*MAT_SPRING_NONLINEAR_ELASTIC',
                   '*MAT_SPRING_GENERAL_NONLINEAR', '*DEFINE_COORDINATE_SYSTEM', '*DEFINE_SD_ORIENTATION',
                   '*CONSTRAINED_EXTRA_NODES_SET', '*SET_NODE_LIST', '*ELEMENT_DISCRETE', '*LOAD_BODY_Z',
                   '*DAMPING_GLOBAL', '*CONTROL_TIMESTEP', '*CONTROL_TERMINATION', '*DATABASE_BINARY_D3PLOT',
		   '*DATABASE_BINARY_D3THDT', '*BOUNDARY_PRESCRIBED_MOTION_NODE', '*BOUNDARY_PRESCRIBED_MOTION_SET'])

  # optional keywords
  optional = set (['*DEFINE_COORDINATE_SYSTEM', '*DEFINE_SD_ORIENTATION', '*CONSTRAINED_EXTRA_NODES_SET',
                   '*SET_NODE_LIST', '*LOAD_BODY_Z', '*DAMPING_GLOBAL', '*BOUNDARY_PRESCRIBED_MOTION_NODE',
		   '*BOUNDARY_PRESCRIBED_MOTION_SET'])

  # keyword field formats: (type, length, optional)
  formats = {'*DEFINE_CURVE': [('int', 10, 0), ('blank', 70, 0), ('int', 10, 1)],
             '*PART': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*PART_INERTIA': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 1),
                               ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                               ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                               ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                               ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                               ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('float', 10, 1)],
             '*MAT_RIGID': [('int', 10, 0), ('int', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                            ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*MAT_SPRING_NONLINEAR_ELASTIC': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*MAT_SPRING_GENERAL_NONLINEAR': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*DEFINE_COORDINATE_SYSTEM': [('int', 10, 0), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                          ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                          ('float', 10, 1), ('float', 10, 1)],
             '*DEFINE_SD_ORIENTATION': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('float', 10, 1),
                                       ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
				       ('float', 10, 1)],
             '*CONSTRAINED_EXTRA_NODES_SET': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*SET_NODE_LIST': [('int', 10, 0), ('blank', 70, 0), ('int', 10, 1)],
             '*ELEMENT_DISCRETE': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0),
                                  ('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*LOAD_BODY_Z': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*DAMPING_GLOBAL': [('int', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                                 ('float', 10, 0), ('float', 10, 0), ('float', 10, 0), ('float', 10, 0),
                                 ('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0)],
             '*CONTROL_TIMESTEP': [('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                  ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                  ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                  ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1)],
             '*CONTROL_TERMINATION': [('float', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                     ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1)],
             '*DATABASE_BINARY_D3PLOT': [('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                        ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                        ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                        ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1)],
             '*DATABASE_BINARY_D3THDT': [('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                        ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1),
                                        ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                        ('int', 10, 1), ('int', 10, 1), ('int', 10, 1), ('int', 10, 1)],
             '*BOUNDARY_PRESCRIBED_MOTION_NODE': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0),
                                                 ('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0),
                                                 ('int', 10, 0), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                                 ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1)],
             '*BOUNDARY_PRESCRIBED_MOTION_SET': [('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0),
                                                ('int', 10, 0), ('int', 10, 0), ('int', 10, 0), ('int', 10, 0),
                                                ('int', 10, 0), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1),
                                                ('float', 10, 1), ('float', 10, 1), ('float', 10, 1), ('float', 10, 1)]}

  # field types
  types = {'int': int, 'float': float, 'blank': str}

  # initialize
  root = {}
  root['PATH'] = path
  root['NODES'] = nodes = {}
  root['SET_NODE_LISTS'] = setnodelists = {}
  line = 0

  # parse
  with open (path, 'r') as f:

    while True:

      l = f.readline (); line += 1

      if not l: break # eof

      l = l.strip()

      if not l: continue # empty line

      if l[0] == '#': continue # comment

      if l[0] != '*': # nodal data
        try:
          n = int(l[0:10])
          x = float(l[10:20])
          y = float(l[20:30])
          z = float(l[30:40])
          nodes[n] = (x, y, z)
        except:
          pass
        continue

      # keyword
      s = l.split(',')
      keyword = s[0].upper()
      
      if keyword not in keywords:
          print('skipping %s' % keyword)
          continue

      # keyword data
      data = root.setdefault (keyword, [])

      # format
      nformatlines = len (formats [keyword])
      format = formats [keyword]

      # read keyword block
      block = []
      nlines = 0
      while nlines < nformatlines:
        l = f.readline (); line += 1

        if not l:
          print('parsing blocks for keyword \"%s\" at line %i failed:' % (keyword, lineno))
          print('%i lines in block but %i lines in format definition - not compatible' % (nlines, nformatlines))
          break # eof

        l = l.rstrip()
        if l and l[0] != '#': # skip comments and empty lines
          block.append (l)
          nlines += 1

      # convert to proper types
      d = {}
      lineno = line - nlines + 1
      for i in range (0, nformatlines):
        dataline = block[i]
        n = len (format[i])
        if n == 3: # (type, length, optional)
          ftype, flen, is_optional = format[i]
        else: # (type, length) - not optional
          ftype, flen = format[i]
          is_optional = 0
        fields = [] # reset for error reporting
        try:
          if ftype == 'blank':
            fieldstr = dataline
            fields.append (types[ftype](fieldstr))
          else:
            j = 0
            while j < len (dataline):
              fieldstr = dataline[j:j+flen]
              j += flen
              if fieldstr and not fieldstr.isspace():
                fields.append (types[ftype](fieldstr))
              elif not is_optional:
                break
        except:
          print('conversion of field \"%s\" in keyword \"%s\" at line %i failed:' % (fname, keyword, lineno))
          print('fields read ok:%s' % fields)
          print('fieldstring:%r, fieldname:%r, type:%r, len:%r, is_optional:%s' % (fieldstr, fname, ftype, flen, is_optional))
          print('remaining line: %r' % dataline)
          raise

        nfields = len (fields)

        if nfields == 0 and is_optional: continue # skip optional empty formats

        for j in range (0, nfields):
          fname = '%s%i' % (ftype, j+1)
          d[fname] = fields[j]

      # special treatment for some keywords
      if keyword == '*SET_NODE_LIST':
        setnodelists [d['int1']] = fields[1:] # first is ID, skip it

      # add to root data
      data.append (d)

  # return root data
  return root

def getcard (root, keyword, **args):
  """ get a specific card from keyword data based on field values """
  
  if keyword in root:
    for card in root[keyword]:
      for field, value in args.items():
        if field in card and card[field] == value: return card

  return None

# assign to keyword data
Keyfile.getcard = getcard