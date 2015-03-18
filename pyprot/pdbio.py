"""
Class inhereted by the `Pdb` base class in `pdbmain`.
Contains methods specialized for PDB file processing.
"""


import urllib
import numpy as np
import pandas as pd
import os

class PdbIO(object):
    def __init__():
        pass
        
    def parse_coordsection(self, dest):
        """
        Parses a PDB file into a pandas DataFrame

        Parameters
        ----------
        dest : `str`.
          Path to the target file. E.g., `"/home/.../desktop/my_pdb.pdb"`
          or list of file contents.
        
        Returns
        ----------
        df : DataFrame.

        """
        rec = [['record', 0, 6],
             ['atomnum', 6, 11],
             ['atomname', 12, 16],
             ['altloc', 17, 18], 
             ['residuename', 17, 20], 
             ['chainid', 21, 22], 
             ['residuenum', 22, 26], 
             ['insertion', 26, 27],
             ['xcoord', 30, 38], 
             ['ycoord', 38, 46], 
             ['zcoord', 47, 54], 
             ['occupancy', 54, 60], 
             ['bfactor', 60, 66],
             ['segmentid', 72, 76], 
             ['element', 76, 78], 
             ['charge', 79, 80]]
             
        coords = []

        if isinstance(dest, str) and os.path.isfile(dest):
            in_file = open(dest, 'r')
        else:
            in_file = dest


        for line in in_file:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            row = []
            for r in rec:
                sline = line.strip()
                try:
                    row.append(sline[r[1]:r[2]].strip())
                except IndexError:
                    pass
            row.append(line)
            coords.append(row)
    
        df = pd.DataFrame(coords, columns=[c[0] for c in rec] + ['origline'])
        df.tail()
    
        if isinstance(dest, str) and os.path.isfile(dest):
            in_file.close()   
    
    
        for c in ('atomnum', 'residuenum', 'segmentid', 'charge'):
            try:
                df[c] = df[c].astype(int)
            except ValueError:
                pass
        for c in ('xcoord', 'ycoord', 'zcoord', 'occupancy', 'bfactor'):
            try:
                df[c] = df[c].astype(float)
            except ValueError:
                pass

        return df

            

    def save_pdb(self, dest):
        """
        Writes the contents of the `Pdb` object stored in `.cont` attribute to PDB file.

        Parameters
        ----------
        dest : `str`.
          Path to the target file. E.g., `"/home/.../desktop/my_pdb.pdb"`
        
        Returns
        ----------

        success : `bool`.
          `True` if file was successfully written.

        """
        try:
            with open(dest, 'w') as out:
                for line in self.cont:
                    out.write(line + '\n')
            success = True
        except IOError as e:
            print(e)
            success = False

        return success

        
    def fetch_rcsb(self, pdb_code):
        """
        Fetches PDB file contents from rcsb.org.

        Parameters
        ----------
        
        pdb_code : `str`.
          A 4-letter PDB code, e.g., `"3eiy"`
        
        Returns
        ----------

        pdb_cont : `list`.
          List of PDB file contents after where list item is a `str` of
          PDB file contents.
            
        """   
        pdb_cont = []
        try:
            response = urllib.request.urlopen('http://www.rcsb.org/pdb/files/%s.pdb' %pdb_code.lower())
            dat = response.read().decode('utf-8')
            pdb_cont = [row.strip() for row in dat.split('\n') if row.strip()]
        except urllib.request.HTTPError as e:
            print('HTTP Error %s' %e.code)
        except urllib.request.URLError as e:
            print('URL Error %s' %e.args) 
        
        return pdb_cont