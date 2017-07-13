def get_include_filename(card_lines, include_dir=''):
  # type: (List[str], str) -> str
  """
  Parses an INCLUDE file split into multiple lines (as a list).

  Parameters
  ----------
  card_lines : List[str]
      the list of lines in the include card (all the lines!)
  include_dir : str; default=''
      the include directory

  Returns
  -------
  filename : str
      the INCLUDE filename
  """
  
  
  
  card_lines2 = []
  
  for line in card_lines:
      line = line.strip('\t\r\n ')
      card_lines2.append(line)

  card_lines2[0] = card_lines2[0][7:].strip()  # truncate the cardname
  filename = ''.join(card_lines2)
  filename = filename.strip('"').strip("'")
  
  
  if ':' in filename:
      ifilepaths = filename.split(':')
      print("           ifilepaths 00    = %s" %(ifilepaths))
      
      filename = os.path.join(*ifilepaths)
      print("           filename_Brut_00    = %s" %(filename))
      print("             LAST PARTIE       = %s" %(ifilepaths[-1]))
  
  
      ## Ajout ACA
      if not os.path.exists(filename) :
        
        print("        Le Fichier %s existe" %(filename) )
        
        
        # Récupération Valeur Variable Environnement Projet
        if nt.environ.has_key(ifilepaths[0].strip()) :
          
          currVarEnv = nt.environ[ifilepaths[0]]
        
          filename = os.path.join(currVarEnv , ifilepaths[1])
          filename_FINAL_02 = currVarEnv + '/' + ifilepaths[1]

          if AcaDebugPrintValeur == 1 :
            print("             currVarEnv     00 = %s" %(currVarEnv))
            print("             filename FINAL 00 = %s" %(filename))
            print("             filename FINAL 02 = %s" %(filename_FINAL_02))
            
            
  
  if is_windows:
      filename = os.path.join(include_dir, filename).replace('/', '\\')
  else:
      filename = os.path.join(include_dir, filename).replace('\\', '/')
  
  


  return filename