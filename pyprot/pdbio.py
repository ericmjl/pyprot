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
        
    def coordsec_to_ary(self, dest):
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
             ['altloc', 16, 17], 
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
             ['charge', 78, 80]]
             
        coords = []

        if isinstance(dest, str) and os.path.isfile(dest):
            in_file = open(dest, 'r')
        else:
            in_file = dest


        for line in in_file:
            if not line.startswith(('ATOM', 'HETATM', 'ANISOU','TER')):
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
                #df[c] = df[c].apply(lambda x: None if not x else x)
                df[c] = df[c].astype(int)
            except ValueError:
                pass
        for c in ('xcoord', 'ycoord', 'zcoord', 'occupancy', 'bfactor'):
            try:
                df[c] = df[c].apply(lambda x: np.nan if not x else x)
                df[c] = df[c].astype(float)
            except ValueError:
                pass

        return df


    def coordsec_to_file(self, dest):
        """
        Writes the contents of the `coord_ary` DataFrame to file.

        Parameters
        ----------
        dest : `str`.
          Path to the target file. E.g., `"/home/.../desktop/my_pdb.pdb"`
        
        """
        df = self.coord_ary.copy()
        df.drop('origline', axis=1, inplace=True)

        df['record'] = df['record'].apply(lambda x: x + (6-len(x))*' ')
        df['atomnum'] = df['atomnum'].apply(lambda x: (5-len(str(x)))*' ' + str(x))
        df['atomname'] = df['atomname'].apply(lambda x: '  ' + str(x) + (3-len(str(x)))*' ' )
        df['residuename'] = df['residuename'].apply(lambda x: str(x) + (3-len(str(x)))*' ' )
        df['altloc'] = df['altloc'].apply(lambda x: ' ' if not x else x)
        df['chainid'] = df['chainid'].apply(lambda x: ' ' + x)
        df['residuenum'] = df['residuenum'].apply(lambda x: ' '*(4-len(str(x))) + str(x))
        df['insertion'] = df['insertion'].apply(lambda x: ' ' if not x else x)
        df['segmentid'] = df['segmentid'].apply(lambda x: ' '*(6-len(str(x))) + str(x))
        df['element'] = df['element'].apply(lambda x: ' '*(6-len(str(x))) + str(x))
        df['charge'] = df['charge'].apply(lambda x: ' '*(2-len(str(x))) + str(x))

        # fix TER entries later
        df_nt = df[df['record'] != 'TER']
        df_nt['xcoord'] = df_nt['xcoord'].apply(lambda x: ' '*(11-len('%.3f' % x)) + '%.3f' % x)
        df_nt['ycoord'] = df_nt['ycoord'].apply(lambda x: ' '*(8-len(str('%.3f' % x))) + '%.3f' % x)
        df_nt['zcoord'] = df_nt['zcoord'].apply(lambda x: ' '*(8-len(str('%.3f' % x))) + '%.3f' % x)
        df_nt['occupancy'] = df_nt['occupancy'].apply(lambda x: ' '*(6-len('%.2f' % x)) + '%.2f' % x)
        df_nt['bfactor'] = df_nt['bfactor'].apply(lambda x: ' '*(6-len('%.2f' % x)) + '%.2f' % x)

        df = df_nt[df_nt.index == df.index]        
        
        # fix TER
        df = df.where((pd.notnull(df)), ' ')
        for r in ('xcoord','ycoord','zcoord','occupancy','bfactor'):
            df.loc[df['record'] == 'TER   ', r] = df.loc[df['record'] == 'TER   ', r]\
                    .apply(lambda x: x.replace('nan', '   '))
                    
        np.savetxt(dest, 
                   df.values, 
                   delimiter='', 
                   newline='\n', 
                   header='', 
                   footer='',
                   fmt='%s',
                   comments='# ')

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