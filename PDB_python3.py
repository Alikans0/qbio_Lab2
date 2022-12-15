#!/usr/bin/env python
# -*- coding: iso-8859-1 -*- 

try:
    import sys
    from math import sqrt, cos, sin, pi, acos
    #from string import rjust
    #from string import ljust
    #import resource
    #sys.path.append("/home/louet/PYTHON-SCRIPT/MOUVEMENTS/")
    #sys.path.append("/home/louet/PYTHON-SCRIPT/")
    #import alignment as ali
    #import matrice as mat
except ImportError:
    print('''\033[31mUn des modules suivants est introuvable:
    sys
    math
    string
    resource\033[0m''')
    sys.exit()

class atome():
    def __init__(self, desc="ATOM", num_atome=1, type_atome="  ", 
                 type_res="   ", chain=" ", num_res=1, x=0.0, y=0.0, z=0.0,
                 occupancy=1.00, bfactor=0.00, segname="         ", 
                 element="   "):
        """Definit un atome tel que decrit dans un PDB
        
        desc = string
        num_atome = integer
        type_atome = string
        type_res = string
        chain = string
        num_res = integer
        x = float
        y = float
        z = float
        occupancy = float
        bfactor = float
        segname = string
        element = string
        """
        self.desc = desc
        self.num_atome = num_atome
        self.type_atome = type_atome
        self.type_res = type_res
        self.chain = chain
        self.num_res = num_res
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.segname = segname
        self.element = element
        self.mass = 1.0
        self.charge = 0.0
        self.radius = 0.0
        self.masses_atomes_adt = {"H": 1.008, "HD": 1.008, "HS": 1.008,
                                  "C": 12.011, "A": 12.011, "N" : 14.007,
                                  "NA": 14.007, "NS": 14.007, "OA": 15.999,
                                  "OS": 15.999, "F": 18.998, "Mg": 24.305,
                                  "MG": 24.305, "P": 30.9738, "SA": 32.06,
                                  "S": 32.06, "Cl": 35.45, "CL": 35.45,
                                  "Ca": 40.08, "CA": 40.08, "Mn": 54.938,
                                  "MN": 54.938, "Fe": 55.847, "FE": 55.847,
                                  "Zn": 65.37, "ZN": 65.37, "Br": 79.904,
                                  "BR": 79.904, "I": 126.90447}
        self.masse_atomes = {"C": 12.011, "O": 15.999, "N": 14.007, "S": 32.060,
                             "H": 1.008, "F": 18.998, "N1+": 14.007}
        self.cgRadii = {'ALA': 1.5, 'ARG': 3.2, 'ASN': 2.0, 'ASP': 2.0,
                        'CYS': 1.8, 'GLN': 2.5, 'GLU': 2.5, 'GLY': 0.0,
                        'HIS': 2.5, 'ILE': 2.0, 'LEU': 2.0, 'LYS': 3.0,
                        'MET': 2.5, 'PHE': 2.7, 'PRO': 1.8, 'SER': 1.8,
                        'THR': 2.8, 'TRP': 3.5, 'TYR': 3.0, 'VAL': 1.8}

    def DisplayAtm(self):
        """Affiche a l'ecran les caracteristique d'un atome
        
        Usage : atome.DisplayAtm()        
        """
        if "ATOM" in self.desc or "HETATM" in self.desc:
            print(">>  SegName: %s\n    Type: %3s\n    Atom index: %5d\n    \
Residu: %3d %s\n    Coordinates: { %.3f  %.3f  %.3f }" %(self.segname,
                                                         self.type_atome,
                                                         self.num_atome,
                                                         self.num_res,
                                                         self.type_res,
                                                         self.x, self.y,
                                                         self.z))

    def DisplayAtmLight(self):
        """Affiche a l'ecran les caracteristiques principales d'un atome
        
        Usage : atome.DisplayAtmLight()        
        """
        print("Type: %4s | Numero de l'atome: %6d | Residu: %4s | \
Numero residu: %4d" %(self.type_atome, self.num_atome,
                      self.type_res, self.num_res))

    def DetermineMass(self):
        try:
            self.mass = self.masse_atomes[self.element]
        except KeyError:
            pass

    def DetermineMassAdt(self):
        try:
            self.mass = self.masse_atomes_adt[self.element]
        except KeyError:
            pass

    def Copy(self):
        """Copie l'atome dans un nouvel obet atom
        
        Usage : NewAtom = atom.Copy()
        """
        new_a = atome(self.desc, self.num_atome, self.type_atome, self.type_res,
                      self.chain, self.num_res, self.x, self.y, self.z,
                      self.occupancy, self.bfactor, self.segname, self.element)
        return new_a


class pdb():
    def __init__(self, fichier="", force=False, verbose=False):
        """Definit un fichier PDB

        file = string
        nom = string
        description = list[description]
        liste_atome = list[atome]
        liste_protein = dico{Code3L : Code1L}
        liste_protein2 = dico{Code1L : Code3L}
        force = Booleen (defaut False - lecture ne s'arrete pas pour erreur)
        verbose = Booleen (defaut True - affiche les avertissements)
        
        """
        self.file = fichier
        self.nom = fichier
        self.description = []
        self.liste_atome = []
        self.force = force
        self.verbose = verbose
        self.connect = ""
        self.liste_protein = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                              'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
                              'HSE': 'H', 'HSD': 'H', 'HIS': 'H', 'ILE': 'I',
                              'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                              'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                              'TYR': 'Y', 'VAL': 'V', 'GLYM': 'G', 'CYSG': 'C',
                              'CYSP': 'C', 'CYG': 'C', 'CYP': 'C', 'GLM': 'G',
                              'GYM': 'G'}
        self.liste_protein2 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
                               'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
                               'H': 'HSE', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU',
                               'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                               'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR',
                               'V': 'VAL', 'H':'HSD', 'G': 'GLYM', 'C': 'CYSG',
                               'C': 'CYSP', 'C': 'CYG', 'C': 'CYP', 'G': 'GLM',
                               'G': 'GYM'}
        self.Main()

    def Main(self):
        """Fonction principale de la classe pdb (lancee au demarrage)"""
        if self.file != "":
            try:
                fichier = open(self.file, 'r')
            except IOError:
                print("Fichier %s introuvable" %(self.file))
                sys.exit()
            else:
                self.StorePDBInfo(fichier.read())

    def MsgErr(self, fonction, message):
        print("\033[31mFonction : " + fonction)
        print(message + "\033[0m")
        print("")
        print("-> Aide " + fonction)
        print("")
        print(getattr(self, fonction).__doc__)
        sys.exit()

    def PDBLoad(self, PDBCode):
        """Charge un PDB depuis rcsb.org

        Usage : pdb.PDBLoad(PDBCode)

        PDBCode = string (4 Lettres)
        """
        try:
            from urllib.request import urlopen
        except ImportError:
            self.MsgErr("PDBLoad", "Le module urllib2 est introuvable !")
        testpdb = isinstance(PDBCode, str)
        if not testpdb or len(PDBCode) != 4:
            self.MsgErr("PDBLoad", "Format argument incorrect !")
        try:
            url = urlopen("https://files.rcsb.org/download/%s.pdb" %PDBCode)
        except:
            self.MsgErr("PDBLoad", "Erreur: Code invalide !")
        self.file = PDBCode
        self.nom = PDBCode
        self.StorePDBInfo(url.read().decode())

    def Copy(self, atom_select):
        """Copie le fichier PDB et le retourne vers un nouveau PDB
        
        Usage : NewPDB = pdb.Copy(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("Copy", "Format argument incorrect !")
        new_pdb = pdb()
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                new_a = atome(a.desc, a.num_atome, a.type_atome, a.type_res,
                              a.chain, a.num_res, a.x, a.y, a.z, a.occupancy,
                              a.bfactor, a.segname, a.element)
                new_pdb.liste_atome.append(new_a)
        return new_pdb

    def DisplayAllAtome(self):
        """Affiche tous les atomes d'un PDB
        
        Usage : pdb.DisplayAllAtome()
        """
        for a in self.liste_atome:
            a.DisplayAtmLight()

    def AjoutAtm(self, atome):
        """Ajoute un atome au PDB
        
        Usage : pdb.AjoutAtm(atome)
        
        atome = atome()
        """
        self.liste_atome.append(atome)

    def DeleteAtm(self, index):
        self.liste_atome.pop(index)

    def AjoutRes(self, code3l, code1l):
        """Ajoute un identifiant de residu a la liste proteique
        
        Usage : pdb.AjoutRes(code3l, code1l)
        
        code3l = string
        code1l = string
        """
        test1l = isinstance(code1l, str)
        test3l = isinstance(code3l, str)
        if not test1l or not test3l:
            self.MsgErr("AjoutRes", "Format argument incorrect !")
        if self.liste_protein.has_key(code3l):
            if self.verbose:
                print("Le r\E9sidu %s existe d\E9j\E0 dans la base de donn\E9e" % code3l)
        else:
            self.liste_protein[code3l] = code1l
            
    def AjoutRes2(self, code3l, code1l):
        """Ajoute un identifiant de residu a la liste proteique
        
        Usage : pdb.AjoutRes(code3l, code1l)
        
        code3l = string
        code1l = string
        """
        test1l = isinstance(code1l, str)
        test3l = isinstance(code3l, str)
        if not test1l or not test3l:
            self.MsgErr("AjoutRes2", "Format argument incorrect !")
        if self.liste_protein2.has_key(code1l):
            if self.verbose:
                print("Le r\E9sidu %s existe d\E9j\E0 dans la base de donn\E9e" % code1l)
        else:
            self.liste_protein[code1l] = code3l

    def CheckSel(self, atom_select, objet):
        """Permet de tester une selection
        
        Usage : pdb.CheckSel(atom_select, objet)
        
        atom_select = string
        objet = string
        """
        selection = self.ParseSelect(atom_select, objet)
        if selection[0]:
            self.MsgErr("ParseSelect", "Selection incorrecte ! (%s)" \
                                    % selection[1])
        else:
            return selection[1]

    def ParseSelect(self, select, objet):
        """Cree une 'chaine de caractere test' en fonction d'une selection
        Sortie = [codeSortie, selection]
            codeSortie = 1 -> erreur ; codeSortie = 0 -> pas d'erreur
            selection = string
        
        Usage : pdb.ParseSelect(select, objet)
        
        select = string
        objet = string
        
        Utilisation de select :
            all      -> tout
            protein  -> residus dont le code3l est stocke dans pdb.liste_protein
            backbone -> = name CA C N O (atomes du backbone)
            name     -> nom atome        (ex: name CA ; name CA C N O HD1)
            index    -> numero atome     (ex: index 1 to 4 7 18 to 30) 
            resid    -> numero residu    (ex: resid 1 to 4 7 18 to 30)
            resname  -> nom du residu    (ex: resname HIS ; resname TRP VAL)
            segname  -> nom du segment   (ex: segname GIA ; segname GIB GIG)
            chain    -> nom de la chaine (ex: chain A ; chain A B)
            occupancy-> nombre           (ex: occupancy 1; occupancy 0)
            
            operateurs logiques : 
                or ; and ; not ; ( ; )
        ATTENTION : Les parentheses doivent comporter des espaces avant et apres
                    Exemple : " ( resid 10 ) and ( name CA ) "
        """
        if len(select) == 0:
            return 1
        split_select = select.split()
        reserve = ["name", "index", "resid", "resname", "segname", "chain",
                   "or", "and", "not", "all", "protein", "(", ")", "water",
                   "backbone"]
        selection_string = "('ATOM' in %s.desc or 'HETATM' in %s.desc or \
'TER' in %s.desc) and " %(objet, objet, objet)
        i = 0
        longueur = len(split_select) - 1
        while i <= longueur:
            if split_select[i] == "index":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve and i <= longueur:
                        if split_select[i + 1] == "to":
                            selection_string += "%s.num_atome >= %s and " \
                                                %(objet, split_select[i])
                            selection_string += "%s.num_atome <= %s or " \
                                                %(objet, split_select[i + 2])
                            i += 3
                        else:
                            selection_string += "%s.num_atome == %s or " \
                                                %(objet, split_select[i])
                            i += 1
                except IndexError:
                    if i != (longueur + 1):
                        selection_string += "%s.num_atome == %s or " \
                                            %(objet, split_select[i])
                        i += 1
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "resid":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve and i <= longueur:
                        if split_select[i + 1] == "to":
                            selection_string += "%s.num_res >= %s and " \
                                                %(objet, split_select[i])
                            selection_string += "%s.num_res <= %s or " \
                                                %(objet, split_select[i + 2])
                            i += 3
                        else:
                            selection_string += "%s.num_res == %s or " \
                                                %(objet, split_select[i])
                            i += 1
                except IndexError:
                    if i != (longueur + 1):
                        selection_string += "%s.num_res == %s or " \
                                            %(objet, split_select[i])
                        i += 1
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "name":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve:
                        selection_string += "'%s' == %s.type_atome or " \
                                            %(split_select[i], objet)
                        i += 1
                except IndexError:
                    pass
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "resname":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve:
                        selection_string += "'%s' in %s.type_res or " \
                                            %(split_select[i], objet)
                        i += 1
                except IndexError:
                    pass
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "segname":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve:
                        selection_string += "'%s' in %s.segname or " \
                                            %(split_select[i], objet)
                        i += 1
                except IndexError:
                    pass
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "chain":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve:
                        selection_string += "%s.chain == '%s' or " \
                                            %(objet, split_select[i])
                        i += 1
                except IndexError:
                    pass
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "occupancy":
                selection_string += "("
                i += 1
                try:
                    while split_select[i] not in reserve:
                        selection_string += "%s.occupancy == %s or " \
                                            %(objet, split_select[i])
                        i += 1
                except IndexError:
                    pass
                if selection_string[-3:-1] == "or":
                    selection_string = selection_string[:-3]
                selection_string += ") "
            elif split_select[i] == "and":
                i += 1
                selection_string += "and "
            elif split_select[i] == "or":
                i += 1
                selection_string += "or "
            elif split_select[i] == "not":
                i += 1
                selection_string += "not "
            elif split_select[i] == "all":
                i += 1
                selection_string += "True "
            elif split_select[i] == "protein":
                i += 1
                selection_string += "("
                for k, v in self.liste_protein.items():
                    selection_string += "%s.type_res == '%s' or " %(objet, k)
                selection_string = selection_string[:-3]
                selection_string += ")"
            elif split_select[i] == "backbone":
                i += 1
                selection_string += "%s.type_atome == 'CA' or %s.type_atome \
== 'C' or %s.type_atome == 'N' or %s.type_atome == 'O'" %(objet, objet, objet,
                                                          objet)
            elif split_select[i] == "(":
                i += 1
                selection_string += "("
            elif split_select[i] == ")":
                i += 1
                selection_string += ") "
            elif split_select[i] == "water":
                i += 1
                selection_string += "%s.type_res == 'HOH' or %s.type_res == \
'TIP3'" %(objet, objet)
            else:
                return [1, select]
        if selection_string[-4:-1] == "and ":
            selection_string = selection_string[:-4]
        if selection_string[-3:-1] == "or":
            selection_string = selection_string[:-3]
        globals()[objet] = atome()
        try:
            eval(selection_string)
        except SyntaxError:
            return [1, select]
        else:
            return [0, selection_string]

    def StorePDBInfo(self, ligne):
        """Lit et stocke les informations contenues dans le fichier PDB
        
        Usage : pdb.StorePDBInfo(ligne)
        
        ligne = string (line = toutes les lignes mises bout a bout du pdb)
        """
        nbLigne = 1
        num_res_hexa = False
        for line in ligne.split("\n"):
            if "ATOM" in line[0:6] or "HETATM" in line[0:6]:
                try:
                    num_atome = int(line[6:11])
                except ValueError:
                    try:
                        num_atome = int("0x" + line[6:11], 16)
                    except ValueError:
                        if self.force:
                            if self.verbose:
                                print("Avertissement: ligne %d (%s) ignoree"\
                                      %(nbLigne, line))
                        else:
                            self.MsgErr("StorePDBInfo", "Erreur lecture %d: %s"\
                                                        %(nbLigne, line))
                try:
                    occupancy = float(line[54:60])
                except ValueError:
                    occupancy = 0.0
                except IndexError:
                    if self.force:
                        if self.verbose:
                            print("Avertissement: ligne %d (%s) ignoree" %(nbLigne,
                                                                           line))
                    else:
                        self.MsgErr("StorePDBInfo", "Erreur lecture %d: %s" \
                                                    %(nbLigne, line))
                try:
                    desc = line[0:6].split()[0]
                    type_atome = line[12:16].split()[0]
                    model = line[16]
                    type_res = line[17:21].split()[0]
                    chain = line[21]
                    #TODO DEALING WITH ANISOU
                    if num_res_hexa:
                        if line[26].isalpha():
                            if line[26] == "A":
                                num_res = int(line[22:26], 16)
                        else:
                            num_res = int(line[22:26], 16)
                    else:
                        if line[26].isalpha():
                            if line[26] == "A":
                                num_res = int(line[22:26])
                        else:
                            num_res = int(line[22:26])
                    x = float(line[27:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    bfactor = float(line[60:66])
                    segname = line[67:76]
                    element = line[76:78]
                except (ValueError, IndexError) as e:
                    if self.force:
                        if self.verbose:
                            print("Avertissement: ligne %d (%s) ignoree" %(nbLigne,
                                                                           line))
                    else:
                        self.MsgErr("StorePDBInfo", "Erreur lecture %d: %s" \
                                                    %(nbLigne, line))
                else:
                    try:
                        a = atome(desc, num_atome, type_atome, type_res, chain, 
                                  num_res, x, y, z, occupancy, bfactor, segname,
                                  element)
                        a.DetermineMass()
                    except UnboundLocalError:
                        print(line)
                        if self.verbose:
                            print("Avertissement: ligne %d (%s) modeles \
multiples non supportes, ignore" %(nbLigne, line))
                    else:
                        if model == "A" or model == " ":
                            self.AjoutAtm(a)
                        else:
                            if self.verbose:
                                print("Avertissement: ligne %d (%s) modeles \
multiples non supportes, ignore" %(nbLigne, line))
            elif "TER" in line[0:3]:
                try:
                    desc = line[0:6].split()[0]
                    num_atome = int(line[6:11])
                    type_res = line[17:21].split()[0]
                    chain = line[21]
                    num_res = int(line[22:26])
                    a = atome(desc, num_atome, "", type_res, chain,
                        num_res)
                except ValueError:
                    if self.force:
                        if self.verbose:
                            print("Avertissement: ligne %d (%s) ignoree" %(nbLigne,
                                                                           line))
                    else:
                        self.MsgErr("StorePDBInfo", "Erreur lecture %d: %s" \
                                                    %(nbLigne, line))
                except IndexError:
                    if self.force:
                        if self.verbose:
                            print("Avertissement: ligne %d (%s) ignoree" %(nbLigne,
                                                                           line))
                    else:
                        self.MsgErr("StorePDBInfo", "Erreur lecture %d: %s" \
                                                    %(nbLigne, line))
                else:
                    self.AjoutAtm(a)
            elif "ENDMDL" in line[0:6]:
                if self.verbose:
                    print("Avertissement: PDB mutli modeles non supportes, \
modele restants ignores")
                break
            elif "CONECT" in line:
                self.connect += line + "\n"
            else:
                self.description.append(line)
            try:
                if num_res >= 9999:
                    num_res_hexa = True
            except UnboundLocalError:
                pass
            nbLigne += 1
        if len(self.liste_atome) == 0:
            print("Erreur: %s aucun atome valide au format PDB" % self.file)
            sys.exit()
        else:
            if self.verbose:
                print("Fichier PDB: %s, lu avec succes" % self.file)

    def AjoutMass(self, psf, atom_select):
        '''Ajoute les masses aux atomes a l'aide d'un fichier psf CHARMM'''
        fipsf = open(psf, 'r')
        stock = False
        i = 0
        for line in fipsf.readlines():
            if "!NBOND" in line:
                break
            if stock and len(line) > 1:
                splitage = line.split()
                self.liste_atome[i].mass = float(splitage[7])
                i += 1
            if "!NATOM" in line:
                stock = True

    def AjoutElement(self, atom_select):
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                try:
                    int(a.type_atome[0])
                except ValueError:
                    a.element = " " + a.type_atome[0]

    def CentreDeMasse(self, atom_select):
        """Calcul le centre de masse d'une selection d'atome(s)
        Sortie = [coordX, coordY, coordZ] (liste)
        
        Usage : pdb.CentreDeMasse(atom_select)
        
        atom_select = string
        """
        bary_X = 0.0
        bary_Y = 0.0
        bary_Z = 0.0
        selection = self.CheckSel(atom_select, "a")
        masseTotale = 0.0
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                try:
                    bary_X += a.x * a.mass
                    bary_Y += a.y * a.mass
                    bary_Z += a.z * a.mass
                    masseTotale += a.mass
                except TypeError:
                    a.DisplayAtm()
                    self.MsgErr("CentreDeMasse", "Format de l'atome incorrect !")
                ok = True
        if not ok:
            self.MsgErr("CentreDeMasse", "Au moins un atome doit etre selectionne !")
        baryX = bary_X / masseTotale
        baryY = bary_Y / masseTotale
        baryZ = bary_Z / masseTotale
        return [baryX, baryY, baryZ]  
        
    def DistanceMatrixR(self, atom_select):
        """Cree la matrice de distance des atomes sous forme de vecteur
        Sortie = matrice (1 ligne), nbatome
        
        Usage : pdb.DistanceMatrixR(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
             self.MsgErr("DistanceMatrixR", "Format argument incorrect !")
        selection = self.CheckSel(atom_select, "a")
        liste_atome = []
        nbatome = 0
        for a in self.liste_atome:
            test = eval(selection)
            if nbatome > 2000:
                self.MsgErr("DistanceMatrixR", "Nbre atome > 2000 pour calcul \
de la matrice !")
            elif test:
                liste_atome.append(a)
                nbatome += 1
        if nbatome <= 1:
            self.MsgErr("DistanceMatrixR", "Nbre atome <= 1 pour calcul de la \
matrice pour le pdb %s" % self.file)
        matrice = []
        for a1 in liste_atome:
            xa = a1.x
            ya = a1.y
            za = a1.z
            for a2 in liste_atome:
                xb = a2.x
                yb = a2.y
                zb = a2.z
                distance = sqrt((xb - xa) * (xb - xa) + (yb - ya) * \
                                (yb - ya) + (zb - za) * (zb - za))
                matrice.append(distance)
        return matrice, len(liste_atome)
        
    def DistanceMatrix(self, atom_select):
        """Cree la matrice de distance des atomes selectionnes
        Sortie = matrice (N x N)
        
        Usage : pdb.DistanceMatrixR(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("DistanceMatrixR", "Format argument incorrect !")
        selection = self.CheckSel(atom_select, "a")
        liste_atome = []
        nbatome = 0
        for a in self.liste_atome:
            test = eval(selection)
            if nbatome > 2000:
                self.MsgErr("DistanceMatrix", "Nbre atome > 2000 pour calcul \
de la matrice !")
            elif test:
                liste_atome.append(a)
                nbatome += 1
        if nbatome <= 1:
            self.MsgErr("DistanceMatrix", "Nbre atome <= 1 pour calcul de la \
matrice !")
        matrice = []
        for a1 in liste_atome:
            xa = a1.x
            ya = a1.y
            za = a1.z
            line = []
            for a2 in liste_atome:
                xb = a2.x
                yb = a2.y
                zb = a2.z
                distance = sqrt((xb - xa) * (xb - xa) + (yb - ya) * \
                                (yb - ya) + (zb - za) * (zb - za))
                line.append(distance)
            matrice.append(line)
        return matrice

    def RenumRes(self, new_num_res, atom_select):
        """Renumerote une selection de residus a partir d'un entier choisi
        
        Usage : pdb.RenumRes(new_num_res, atom_select)
        
        new_num_res = integer
        atom_select = string
        """
        testnum = isinstance(new_num_res, int)
        testsel = isinstance(atom_select, str)
        if not testnum or not testsel:
            self.MsgErr("RenumRes", "Format argument(s) incorrect(s) !")
        new_num_res -= 1
        current_res = -1
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                if current_res != a.num_res:
                    current_res = a.num_res
                    new_num_res += 1
                    a.num_res = new_num_res
                else:
                    a.num_res = new_num_res

    def RenumAtm(self, new_num_atom, atom_select):
        """Renumerote une selection d'atomes a partir d'un entier choisi
        
        Usage : pdb.RenumAtm(new_num_atom, atom_select)
        
        new_num_atom = integer
        atom_select = string
        """
        testnum = isinstance(new_num_atom, int)
        testsel = isinstance(atom_select, str)
        if not testnum or not testsel:
            self.MsgErr("RenumAtm", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.num_atome = new_num_atom
                new_num_atom += 1

    def GetResNameSeq(self, atom_select):
        res_actuel = -1
        seq = []
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                if res_actuel != a.num_res:
                    res_actuel = a.num_res
                    seq.append(a.type_res)
        return seq
        
    def GetSeq(self, atom_select):
        res_actuel = -1
        seq = ""
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                if res_actuel != a.num_res:
                    res_actuel = a.num_res
                    seq += self.liste_protein[a.type_res]
        return seq
        
    def GetResid(self, atom_select):
        res_actuel = -1
        listResid = []
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                if res_actuel != a.num_res:
                    res_actuel = a.num_res
                    listResid.append(a.num_res)
        return listResid

    def GetCoorSeq(self, atom_select):
        seq = []
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                seq.append([a.x, a.y, a.z])
        return seq

    def ModifChaine(self, chain, atom_select):
        """Modifie la chaine de la selection d'atome
        
        Usage : pdb.ModifChaine(chain, atom_select)
        
        chain = string
        atom_select = string
        """
        testchain = isinstance(chain, str)
        testsel = isinstance(atom_select, str)
        if not testchain or not testsel:
            self.MsgErr("ModifChaine", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            try:
                test = eval(selection)
                if test:
                    a.chain = chain
            except TypeError:
                pass

    def ModifSegName(self, segname, atom_select):
        """Modifie le nom de segment de la selection d'atome
        
        Usage : pdb.ModifSegName(segname, atom_select)
        
        segname = string
        atom_select = string
        """
        testsegname = isinstance(segname, str)
        testsel = isinstance(atom_select, str)
        if not testsegname or not testsel:
            self.MsgErr("ModifSegName", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.segname = segname

    def ModifBfac(self, bfactor, atom_select):
        """Modifie le B-Factor de la selection d'atome
        
        Usage : pdb.ModifBfac(bfactor, atom_select)
        
        bfactor = float
        atom_select = string
        """
        testbfac = isinstance(bfactor, int) or isinstance(bfactor, float)
        testsel = isinstance(atom_select, str)
        if not testbfac or not testsel:
            self.MsgErr("ModifBfac", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.bfactor = bfactor

    def ModifOccupancy(self, occupancy, atom_select):
        """Modifie l'occupancy de la selection d'atome
        
        Usage : pdb.ModifOccupancy(occupancy, atom_select)
        
        occupancy = float
        atom_select = string
        """
        testocc = isinstance(occupancy, int) or isinstance(occupancy, float)
        testsel = isinstance(atom_select, str)
        if not testocc or not testsel:
            self.MsgErr("ModifOccupancy",
                        "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.occupancy = occupancy

    def ModifResName(self, resname, atom_select):
        """Modifie le nom de residu de la selection d'atome
        
        Usage : pdb.ModifResName(resname, atom_select)
        
        resname = string
        atom_select = string
        """
        testresname = isinstance(resname, str)
        testsel = isinstance(atom_select, str)
        if not testresname or not testsel:
            self.MsgErr("ModifResName", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.type_res = resname

    def ModifAtmName(self, atmname, atom_select):
        """Modifie le nom d'atome de la selection d'atome
        
        Usage : pdb.ModifAtmName(atmname, atom_select)
        
        atmname = string
        atom_select = string
        """
        testatm = isinstance(atmname, str)
        testsel = isinstance(atom_select, str)
        if not testatm or not testsel:
            self.MsgErr("ModifAtmName", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.type_atome = atmname

    def ModifDesc(self, desc, atom_select):
        """Modifie la description de la selection d'atome (ATOM ; HETATM ; TER)
        
        Usage : pdb.ModifDesc(desc, atom_select)
        
        desc = string
        atom_select = string
        """
        testdesc = isinstance(desc, str)
        testsel = isinstance(atom_select, str)
        if not testdesc or not testsel:
            self.MsgErr("ModifDesc", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.desc = desc

    def ModifElem(self, elem, atom_select):
        """Modifie le nom de l'element de la selection d'atome
        
        Usage : pdb.ModifElem(elem, atom_select)
        
        elem = string
        atom_select = string
        """
        testelem = isinstance(elem, str)
        testsel = isinstance(atom_select, str)
        if not testelem or not testsel:
            self.MsgErr("ModifElem", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.element = elem

    def ModifSequence(self, sequence):
        """Remplace la sequence du PDB en fonction d'une autre sequence
        
        Usage : pdb.ModifChaine(sequence)
        
        sequence = string (code1L sans espace ni retour a la ligne)
        """
        for a in self.liste_atome:
            try:
                a.resname = self.liste_protein2[sequence[a.num_res - 1]]
            except KeyError:
                self.MsgErr("ModifSequence", "Residu %s inconnu !\nPour l'ajouter\
, Utiliser pdb.AjoutRes2(code3l, code1l)" % sequence[a.num_res - 1])

    def Translation(self, x, y, z, atom_select):
        """Translate la selection d'atome selon les 3 coordonnees cartesiennes
        
        Usage : pdb.Translation(x, y, z, atom_select)
        
        x = y = z = float
        atom_select = string
        """
        testx = isinstance(x, float) or isinstance(x, int)
        testy = isinstance(y, float) or isinstance(y, int)
        testz = isinstance(z, float) or isinstance(z, int)
        testsel = isinstance(atom_select, str)
        if not testx or not testy or not testz or not testsel:
            self.MsgErr("Translation", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.x += x
                a.y += y
                a.z += z

    def Rotation(self, x, y, z, atom_select):
        """Rote la selection d'atome selon les 3 axes cartesiens
        
        Usage : pdb.Rotation(x, y, z, atom_select)
        
        x = y = z = float (angles en degre)
        atom_select = string
        """
        testx = isinstance(x, float) or isinstance(x, int)
        testy = isinstance(y, float) or isinstance(y, int)
        testz = isinstance(z, float) or isinstance(z, int)
        testsel = isinstance(atom_select, str)
        if not testx or not testy or not testz or not testsel:
            self.MsgErr("Rotation", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                # Rotation selon X
                if x != 0.0:
                    angle = (pi * x) / 180
                    coory = a.y
                    coorz = a.z
                    a.y = (coory * cos(angle)) - (coorz * sin(angle))
                    a.z = (coory * sin(angle)) + (coorz * cos(angle))
                # Rotation selon Y
                if y != 0.0:
                    angle = (pi * y) / 180
                    coorx = a.x
                    coorz = a.z
                    a.x = (coorx * cos(angle)) + (coorz * sin(angle))
                    a.z = -(coorx * sin(angle)) + (coorz * cos(angle))
                # Rotation selon Z
                if z != 0.0:
                    angle = (pi * z) / 180
                    coorx = a.x
                    coory = a.y
                    a.x = (coorx * cos(angle)) - (coory * sin(angle))
                    a.y = (coorx * sin(angle)) + (coory * cos(angle))

    def Center(self, atom_select):
        """Le barycentre de la selection d'atome est centre en 0 (bouge tout)
        
        Usage : pdb.Center(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("Center", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        barycentre = self.CentreDeMasse(atom_select)
        baryX = barycentre[0]
        baryY = barycentre[1]
        baryZ = barycentre[2]
        for a in self.liste_atome:
            a.x -= baryX
            a.y -= baryY
            a.z -= baryZ

    def OrientationZ(self, atom_select):
        """Oriente en Z la selection d'atome (>2) et la centre en 0 (bouge tout)
        
        Usage : pdb.OrientationZ(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("OrientationZ", "Format argument(s) incorrect(s) !")
        matrice = self.DistanceMatrix(atom_select)
        maximum = matrice[0][0]
        coordmax = [0, 0]
        i = 0
        for line in matrice:
            j = 0
            for valeur in line:
                if valeur > maximum:
                    maximum = valeur
                    coordmax = [i, j]
                j += 1
            i += 1
        liste_test = []
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                liste_test.append(a)
        atome1 = liste_test[coordmax[0]]
        atome2 = liste_test[coordmax[1]]
        xa = atome1.x
        ya = atome1.y
        za = atome1.z
        xb = atome2.x
        yb = atome2.y
        zb = atome2.z
        #TRANSLATION
        if xa > xb:
            x = xa - abs(xb - xa) / 2
        else:
            x = xb - abs(xb - xa) / 2
        if ya > yb:
            y = ya - abs(yb - ya) / 2
        else:
            y = yb - abs(yb - ya) / 2
        if za > zb:
            z = za - abs(zb - za) / 2
        else:
            z = zb - abs(zb - za) / 2
        self.Translation(-x, -y, -z, "all")
        #ROTATION
        iteration = 1
        testxa = abs(atome1.x)
        testya = abs(atome2.y)
        while (testxa >= 0.1 or testya >= 0.1) and iteration <= 10:
            xa = atome1.x
            ya = atome1.y
            za = atome1.z
            distance_xyz = sqrt((xa * xa) + (ya * ya) + (za * za))
            angle_rad_y = acos(ya / distance_xyz)
            angle_deg_y = 180 * (angle_rad_y) / pi
            self.Rotation(0.0, 0.0, angle_deg_y, "all")
            za = atome1.z
            angle_rad_z = acos(za / distance_xyz)
            angle_deg_z = 180 * (angle_rad_z) / pi
            self.Rotation(angle_deg_z, 0.0, 0.0, "all")
            iteration += 1
            testxa = abs(atome1.x)
            testya = abs(atome2.y)
        self.Center(atom_select)

    def OrientationGrid(self, atom_select):
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("OrientationGrid", "Format argument(s) incorrect(s) !")
        matrice = self.DistanceMatrix(atom_select)
        liste_atome = []
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                liste_atome.append(a)
        self.Center(atom_select)
        self.ExportPDB("centered.pdb", "all")  
        ### Matrice en X
        matriceX = []
        for a1 in liste_atome:
            ligne = []
            xa = a1.x
            ligne = []
            for a2 in liste_atome:
                xb = a2.x
                distance = sqrt((xb - xa) * (xb - xa))
                ligne.append(distance)
            matriceX.append(ligne)
        maxX = matriceX[0][0]
        coordmaxX = [0, 0]
        i = 0
        for line in matriceX:
            j = 0
            for valeur in line:
                if valeur > maxX:
                    maxX = valeur
                    coordmaxX = [i, j]
                j += 1
            i += 1
        print("X", liste_atome[coordmaxX[0]].num_atome, liste_atome[coordmaxX[1]].num_atome)
        atome1X = liste_atome[coordmaxX[0]]
        atome2X = liste_atome[coordmaxX[1]]
        xa = atome1X.x
        ya = atome1X.y
        za = atome1X.z
        xb = atome2X.x
        yb = atome2X.y
        zb = atome2X.z
        iteration = 1
        testya = abs(atome1X.y)
        testza = abs(atome2X.z)
        while (testya >= 0.1 or testza >= 0.1) and iteration <= 10:
            xa = atome1X.x
            ya = atome1X.y
            za = atome1X.z
            distance_xyz = sqrt((xa * xa) + (ya * ya) + (za * za))
            angle_rad_z = acos(za / distance_xyz)
            angle_deg_z = 180 * (angle_rad_z) / pi
            self.Rotation(angle_deg_z, 0.0, 0.0, "all")
            xa = atome1X.x
            angle_rad_x = acos(xa / distance_xyz)
            angle_deg_x = 180 * (angle_rad_x) / pi
            self.Rotation(0.0, angle_deg_x, 0.0, "all")
            iteration += 1
            testxa = abs(atome1X.x)
            testya = abs(atome2X.y)
        self.Center(atom_select)
        self.ExportPDB("rotationX.pdb", "all")
        ### Matrice en Y
        matriceY = []
        for a1 in liste_atome:
            ligne = []
            ya = a1.y
            ligne = []
            for a2 in liste_atome:
                yb = a2.y
                distance = sqrt((yb - ya) * (yb - ya))
                ligne.append(distance)
            matriceY.append(ligne)
        maxY = matriceY[0][0]
        coordmaxY = [0, 0]
        i = 0
        for line in matriceY:
            j = 0
            for valeur in line:
                if valeur > maxY:
                    maxY = valeur
                    coordmaxY = [i, j]
                j += 1
            i += 1
        print("Y", liste_atome[coordmaxY[0]].num_atome, liste_atome[coordmaxY[1]].num_atome)
        atome1Y = liste_atome[coordmaxY[0]]
        atome2Y = liste_atome[coordmaxY[1]]
        xa = atome1Y.x
        ya = atome1Y.y
        za = atome1Y.z
        xb = atome2Y.x
        yb = atome2Y.y
        zb = atome2Y.z
        iteration = 1
        testxa = abs(atome1Y.x)
        testza = abs(atome2Y.z)
        while (testxa >= 0.1 or testza >= 0.1) and iteration <= 10:
            xa = atome1Y.x
            ya = atome1Y.y
            za = atome1Y.z
            distance_xyz = sqrt((xa * xa) + (ya * ya) + (za * za))
            #angle_rad_x = acos(xa / distance_xyz)
            #angle_deg_x = 180 * (angle_rad_x) / pi
            #self.Rotation(0.0, angle_deg_x, 0.0, "all")
            #ya = atome1Y.y
            angle_rad_y = acos(ya / distance_xyz)
            angle_deg_y = 180 * (angle_rad_y) / pi
            if atome1Y.x > 0:
                self.Rotation(angle_deg_y * -1, 0.0, 0.0, "all")
            else:
                self.Rotation(angle_deg_y, 0.0, 0.0, "all")
            iteration += 1
            testxa = abs(atome1Y.x)
            testya = abs(atome2Y.y)
        self.Center(atom_select)
        self.ExportPDB("rotationY.pdb", "all")
        """
        ### Matrice en Z
        matriceZ = []
        for a1 in liste_atome:
            ligne = []
            za = a1.z
            ligne = []
            for a2 in liste_atome:
                zb = a2.z
                distance = sqrt((zb - za) * (zb - za))
                ligne.append(distance)
            matriceZ.append(ligne)
        maxZ = matriceZ[0][0]
        coordmaxZ = [0, 0]
        i = 0
        for line in matriceZ:
            j = 0
            for valeur in line:
                if valeur > maxZ:
                    maxZ = valeur
                    coordmaxZ = [i, j]
                j += 1
            i += 1
        print "Z", liste_atome[coordmaxZ[0]].num_atome, liste_atome[coordmaxZ[1]].num_atome
        atome1Z = liste_atome[coordmaxZ[0]]
        atome2Z = liste_atome[coordmaxZ[1]]
        xa = atome1Z.x
        ya = atome1Z.y
        za = atome1Z.z
        xb = atome2Z.x
        yb = atome2Z.y
        zb = atome2Z.z
        iteration = 1
        testxa = abs(atome1Z.x)
        testya = abs(atome2Z.y)
        while (testxa >= 0.1 or testya >= 0.1) and iteration <= 10:
            xa = atome1Z.x
            ya = atome1Z.y
            za = atome1Z.z
            distance_xyz = sqrt((xa * xa) + (ya * ya) + (za * za))
            angle_rad_y = acos(ya / distance_xyz)
            angle_deg_y = 180 * (angle_rad_y) / pi
            self.Rotation(0.0, 0.0, angle_deg_y, "all")
            za = atome1Z.z
            angle_rad_z = acos(za / distance_xyz)
            angle_deg_z = 180 * (angle_rad_z) / pi
            self.Rotation(angle_deg_z, 0.0, 0.0, "all")
            iteration += 1
            testxa = abs(atome1Z.x)
            testya = abs(atome2Z.y)
        self.Center(atom_select)
        self.ExportPDB("rotationZ.pdb", "all")
        """

    def Vel2Gromacs(self, atom_select):
        """Divise par 10 les colonnes de x y et z pour un fichier de vitesse
        pour qu'il soit compatible avec gromacs
        
        Usage : pdb.Vel2Gromacs(atom_select)
        
        atom_select = string
        """
        testsel = isinstance(atom_select, str)
        if not testsel:
            self.MsgErr("Vek2Gromacs", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                a.x = a.x / 10
                a.y = a.y / 10
                a.z = a.z / 10

    def ExportFasta(self, fichier_out, atom_select):
        """Exporte au format fasta
        
        Usage : pdb.ExportFasta(fichier_out, atom_select)
        
        fichier_out = string
        atom_select = string
        """
        testout = isinstance(fichier_out, str)
        testsel = isinstance(atom_select, str)
        if not testout or not testsel:
            self.MsgErr("ExportFasta", "Format argument(s) incorrect(s) !")
        fichier = open(fichier_out, 'w')
        ligne = 80
        pas = 1
        res_actuel = -1
        selection = self.CheckSel(atom_select, "a")
        sequence = ""
        nbRes = 0
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                if res_actuel != a.num_res:
                    res_actuel = a.num_res
                    try:
                        sequence += self.liste_protein[a.type_res]
                    except KeyError:
                        self.MsgErr("ExportFasta", "Residu %s inconnu !\n \
Pour l'ajouter, utiliser pdb.AjoutRes(code3l, code1l)" % a.type_res)
                    else:
                        nbRes += 1
                    if (pas // ligne) == 1:
                        ligne += 80
                        sequence += "\n"
                    pas += 1
        enTete = "> src=%s:nbres=%d\n" %(self.nom, nbRes)
        fichier.write(enTete + sequence + "\n")
        fichier.close()
        if self.verbose:
            print("Fichier FASTA sauvegarde dans: %s" % fichier_out)

    def ExportXYZ(self, fichier_out, atom_select):
        """Exporte au format XYZ
        
        Usage : pdb.ExportXYZ(fichier_out, atom_select)
        
        fichier_out = string
        atom_select = string
        """
        testout = isinstance(fichier_out, str)
        testsel = isinstance(atom_select, str)
        if not testout or not testsel:
            self.MsgErr("ExportXYZ", "Format argument(s) incorrect(s) !")
        selection = self.CheckSel(atom_select, "a")
        fichier = open(fichier_out, 'w')
        ligne = ""
        nbAtm = 0
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                x = str(a.x).rjust(15)
                y = str(a.y).rjust(15)
                z = str(a.z).rjust(15)
                ligne += "%3s    %s%s%s\n" %(a.type_atome, x, y, z)
                nbAtm += 1
        fichier.write("%d\n%s -> %s\n" %(nbAtm, self.file, fichier_out) + ligne)
        fichier.close()

    def ExportPDB(self, fichier_out, atom_select="all", TER=False, large=False,
                  connect=False, renumAtm=True):
        """Exporte au format PDB
        
        Usage : pdb.ExportPDB(fichier_out, atom_select)
        
        fichier_out = string
        atom_select = string
        TER = booleen (ecrit le tag TER a chaque fin de segment) -> defaut False
        """
        testout = isinstance(fichier_out, str)
        testsel = isinstance(atom_select, str)
        if not testout or not testsel:
            self.MsgErr("ExportPDB", "Format argument(s) incorrect(s) !")
        if len(self.liste_atome) < 1:
            self.MsgErr("ExportPDB", "Le PDB ne contient aucun atome !")
        selection = self.CheckSel(atom_select, "a")
        fichier = open(fichier_out, 'w')
        nbAtm = 1
        a = self.liste_atome[0]
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                addTER = False
                if renumAtm:
                    if nbAtm >= 100000:
                        num_atome = hex(nbAtm)[2:]
                    else:
                        num_atome = str(nbAtm)
                else:
                    if a.num_atome >= 100000:
                        num_atome = hex(a.num_atome)[2:]
                    else:
                        num_atome = str(a.num_atome)
                if ((nbAtm > 1 and (a.segname != last.segname or \
                   a.chain != last.chain) and TER) or\
                   ("TER" in a.desc and TER)) and desclast == "ATOM":
                    type_res = last.type_res.ljust(4)
                    ligne = "%-6s%5s      %4s%1s%4d\n" %("TER", num_atome,
                                                         type_res, last.chain,
                                                         last.num_res)
                    addTER = True
                    fichier.write(ligne)
                    nbAtm += 1
                    if renumAtm:
                        if nbAtm >= 100000:
                            num_atome = hex(nbAtm)[2:]
                        else:
                            num_atome = str(nbAtm)
                    else:
                        if nbAtm >= 100000:
                            num_atome = hex(a.num_atome)[2:]
                        else:
                            num_atome = str(a.num_atome)
                if "ATOM" in a.desc or "HETATM" in a.desc:
                    if a.bfactor <= -100000 or a.bfactor >= 1000000:
                        if self.verbose:
                            print("Attention BFactor trop grand... ignore...")
                        bfac = "0.00".rjust(6)
                    if a.bfactor <= -1000 or a.bfactor >= 10000:
                        bfac = ("%d" %(int(a.bfactor))).rjust(6)
                    if a.bfactor <= -100 or a.bfactor >= 1000:
                        bfac = ("%3.1f" %(a.bfactor)).rjust(6)
                    else:
                        bfac = ("%3.2f" %(a.bfactor)).rjust(6)
                    x = ("%4.3f" %(a.x)).rjust(8)
                    y = ("%4.3f" %(a.y)).rjust(8)
                    z = ("%4.3f" %(a.z)).rjust(8)
                    type_res = a.type_res.ljust(4)
                    if a.num_res >= 10000:
                        num_res = hex(a.num_res)[2:]
                    else:
                        num_res = str(a.num_res)
                    if len(a.type_atome) == 4:
                        type_atome = a.type_atome.ljust(5)
                        ligne = "%-6s%5s %4s%4s%1s%4s    %-s%-s%-s  %3.2f%-6s\
 %9s%s\n" %(a.desc, num_atome, type_atome, type_res, a.chain, num_res, x,
                y, z, a.occupancy, bfac, a.segname, a.element)
                        fichier.write(ligne)
                        nbAtm += 1
                    else:
                        type_atome = a.type_atome.ljust(4)
                        #MODIFICATION Large atom number
                        if large:
                            num_atome = num_atome.ljust(6)
                            ligne = "%-6s%6s %4s%4s%1s%4s    %-s%-s%-s  %3.2f\
%-6s %9s%s\n" %(a.desc, num_atome, type_atome, type_res, a.chain,
                    num_res, x, y, z, a.occupancy, bfac, a.segname, a.element)
                        #FIN MODIFICATION
                        else:
                            ligne = "%-6s%5s  %4s%4s%1s%4s    %-s%-s%-s  %3.2f\
%-6s %9s%s\n" %(a.desc, num_atome, type_atome, type_res, a.chain,
                    num_res, x, y, z, a.occupancy, bfac, a.segname, a.element)
                        fichier.write(ligne)
                        nbAtm += 1
                if addTER:
                    desclast = "TER"
                else:
                    desclast = a.desc
                last = a
        if not "TER" in a.desc and TER and a.desc == "ATOM" and not addTER:
            if nbAtm >= 100000:
                num_atome = hex(nbAtm)[2:]
            else:
                num_atome = str(nbAtm)
            if last.num_res >= 10000:
                num_res = hex(last.num_res)[2:]
            else:
                num_res = str(last.num_res)
            type_res = last.type_res.ljust(4)
            ligne = "%-6s%5s      %4s%1s%4s\n" %("TER", num_atome, type_res,
                                                 last.chain, num_res)
            fichier.write(ligne)
        if connect and self.connect != "":
            fichier.write(self.connect)
        fichier.write("END\n")
        fichier.close()
        if self.verbose:
            print("Fichier PDB sauvegarde dans: %s" % fichier_out)

    def ExportPQR(self, fichier_out, atom_select="all", TER=False):
        """Exporte au format PQR
        
        Usage : pdb.ExportPQR(fichier_out, atom_select)
        
        fichier_out = string
        atom_select = string
        TER = booleen (ecrit le tag TER a chaque fin de segment) -> defaut False
        """
        testout = isinstance(fichier_out, str)
        testsel = isinstance(atom_select, str)
        if not testout or not testsel:
            self.MsgErr("ExportPDB", "Format argument(s) incorrect(s) !")
        if self.liste_atome < 1:
            self.MsgErr("ExportPDB", "Le PDB ne contient aucun atome !")
        selection = self.CheckSel(atom_select, "a")
        fichier = open(fichier_out, 'w')
        nbAtm = 1
        a = self.liste_atome[0]
        for a in self.liste_atome:
            test = eval(selection)
            if test:
                addTER = False
                if nbAtm >= 100000:
                    num_atome = hex(nbAtm)[2:]
                else:
                    num_atome = str(nbAtm)
                if ((nbAtm > 1 and (a.segname != last.segname or \
                   a.chain != last.chain) and TER) or\
                   ("TER" in a.desc and TER)) and desclast == "ATOM":
                    type_res = last.type_res.ljust(4)
                    ligne = "%-6s%5s      %4s%1s%4d\n" %("TER", num_atome,
                                                         type_res, last.chain,
                                                         last.num_res)
                    addTER = True
                    fichier.write(ligne)
                    nbAtm += 1
                    if nbAtm >= 100000:
                        num_atome = hex(nbAtm)[2:]
                    else:
                        num_atome = str(nbAtm)
                if "ATOM" in a.desc or "HETATM" in a.desc:
                    if a.charge <= -1000000 or a.charge >= 10000000:
                        if self.verbose:
                            print("Attention charge trop grande... ignore...")
                        charge = "0.00".rjust(7)
                    if a.charge <= -10000 or a.charge >= 100000:
                        charge = ("%d" %(int(a.charge))).rjust(7)
                    if a.charge <= -1000 or a.charge >= 10000:
                        charge = ("%3.2f" %(a.charge)).rjust(7)
                    if a.charge <= -100 or a.charge >= 1000:
                        charge = ("%3.3f" %(a.charge)).rjust(7)
                    else:
                        charge = ("%3.4f" %(a.charge)).rjust(7)
                    if a.radius <= -1000000 or a.radius >= 10000000:
                        if self.verbose:
                            print("Attention radius trop grand... ignore...")
                        radius = "0.00".rjust(7)
                    if a.radius <= -10000 or a.radius >= 100000:
                        radius = ("%d" %(int(a.radius))).rjust(7)
                    if a.radius <= -1000 or a.radius >= 10000:
                        radius = ("%3.2f" %(a.radius)).rjust(7)
                    if a.radius <= -100 or a.radius >= 1000:
                        radius = ("%3.3f" %(a.radius)).rjust(7)
                    else:
                        radius = ("%3.4f" %(a.radius)).rjust(7)
                    x = ("%4.3f" %(a.x)).rjust(8)
                    y = ("%4.3f" %(a.y)).rjust(8)
                    z = ("%4.3f" %(a.z)).rjust(8)
                    type_res = a.type_res.ljust(4)
                    if a.num_res >= 10000:
                        num_res = hex(a.num_res)[2:]
                    else:
                        num_res = a.num_res
                    if len(a.type_atome) == 4:
                        type_atome = a.type_atome.ljust(5)
                        ligne = "%-6s%5s %4s%4s%1s%4s    %-s%-s%-s %-7s%-7s\
\n" %(a.desc, num_atome, type_atome, type_res, a.chain, num_res, x, y, z,
      charge, radius)
                        fichier.write(ligne)
                        nbAtm += 1
                    else:
                        type_atome = a.type_atome.ljust(4)
                        ligne = "%-6s%5s  %4s%4s%1s%4s    %-s%-s%-s %-7s%-7s\
\n" %(a.desc, num_atome, type_atome, type_res, a.chain, num_res, x, y, z,
      charge, radius)
                        fichier.write(ligne)
                        nbAtm += 1
                if addTER:
                    desclast = "TER"
                else:
                    desclast = a.desc
                last = a
        if not "TER" in a.desc and TER and a.desc == "ATOM" and not addTER:
            if nbAtm >= 100000:
                num_atome = hex(nbAtm)[2:]
            else:
                num_atome = str(nbAtm)
            if last.num_res > 10000:
                num_res = hex(last.num_res)[2:]
            else:
                num_res = last.num_res
            type_res = last.type_res.ljust(4)
            ligne = "%-6s%5s      %4s%1s%4s\n" %("TER", num_atome, type_res,
                                                 last.chain, num_res)
            fichier.write(ligne)
        fichier.write("END\n")
        fichier.close()
        if self.verbose:
            print("Fichier PQR sauvegarde dans: %s" % fichier_out)

    def AddChargeRadiusFromPqr(self, pqrFile, atom_select):
        fi = open(pqrFile, 'r')
        dico = {}
        for line in fi.readlines():
            if "ATOM" in line[0:6]:
                type_atome = line[12:16].split()[0]
                type_res = line[17:21].split()[0]
                charge = float(line[54:62])
                radius = float(line[62:68])
                dico[type_atome + "_" + type_res] = [charge, radius]
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            try:
                a.charge = dico[a.type_atome + "_" + a.type_res][0]
                a.radius = dico[a.type_atome + "_" + a.type_res][1]
            except KeyError:
                pass

    def AddChargeRadiusFromPrmFile(self, prmFile, atom_select):
        """Prm file structured as follow:
        #Residue_name Atom_name Charge Radius
        """
        fi = open(prmFile, 'r')
        dico = {}
        for line in fi.readlines():
            if line[0] != '#':
                l = line.split()
                dico[l[1] + "_" + l[0]] = [float(l[2]), float(l[3])]
        selection = self.CheckSel(atom_select, "a")
        for a in self.liste_atome:
            test = eval(selection)
            try:
                a.charge = dico[a.type_atome + "_" + a.type_res][0]
                a.radius = dico[a.type_atome + "_" + a.type_res][1]
            except KeyError:
                pass
  
##############################################################################
###### MATRICE DEPLACEMENT ###################################################
##############################################################################
def DistanceMatrix(PDB, atom_select):
    '''Retourne la matrice de distance d'une selection d'atome'''
    selection = PDB.CheckSel(atom_select, "a")
    liste_atome = []
    for a in PDB.liste_atome:
        test = eval(selection)
        if test:
            liste_atome.append(a)
    matrice = []
    for a1 in liste_atome:
        ligne = []
        xa = a1.x
        ya = a1.y
        za = a1.z
        ligne = []
        for a2 in liste_atome:
            xb = a2.x
            yb = a2.y
            zb = a2.z
            distance = sqrt((xb - xa) * (xb - xa) + (yb - ya) * \
                            (yb - ya) + (zb - za) * (zb - za))
            ligne.append(distance)
        matrice.append(ligne)
    return matrice


def ExportMatriceDeplacement(PDB1, PDB2, atom_select, out):
    '''Export dans un fichier de sortie la matrice de deplacement de deux PDB'''
    try:
        fichier_out = open(out, 'w')
    except IOError:
        print("\033[31mErreur: Impossible d'ecrire le fichier %s\033[0m"\
              %(out))
        sys.exit()
    matrice1 = DistanceMatrix(PDB1, atom_select)
    matrice2 = DistanceMatrix(PDB2, atom_select)
    for i in range(len(matrice1)):
        for j in range(len(matrice1)):
            ecrire = "%.3f" %(matrice1[i][j] - matrice2[i][j])
            fichier_out.write("%9s"%(ecrire))
        fichier_out.write("\n")   
    fichier_out.close()


def MatriceDeplacement(PDB1, PDB2, atom_select):
    '''Cree une matrice de deplacement de deux PDB'''
    matrice1 = DistanceMatrix(PDB1, atom_select)
    matrice2 = DistanceMatrix(PDB2, atom_select)
    if len(matrice1) != len(matrice2):
        print("\033[31mErreur: Les PDB ne comportent %s pas le meme nombre \
d'atomes apres selection\033[0m")
    matriceres = []
    for i in range(len(matrice1)):
        for j in range(len(matrice1)):
            matriceres.append(matrice1[i][j] - matrice2[i][j])
    return matriceres


def CoeffCorrelationCompletR(PDB1, PDB2, PDB3, PDB4, atom_select):
    '''Calcul le coefficient de correlation entre deux mouvements (-> 4 PDB)'''
    try:
        from rpy2.robjects.packages import importr
        import rpy2.robjects as robjects
    except ImportError:
        print("\033[31mFonction : CoeffCorrelationCompletR")
        print("Module rpy2 requis !\033[0m")
        sys.exit()
    (matrice1, nb_atome1) = PDB1.DistanceMatrixR(atom_select)
    (matrice2, nb_atome2) = PDB2.DistanceMatrixR(atom_select)
    (matrice3, nb_atome3) = PDB3.DistanceMatrixR(atom_select)
    (matrice4, nb_atome4) = PDB4.DistanceMatrixR(atom_select)
    if nb_atome1 != nb_atome2 or nb_atome1 != nb_atome3 or \
       nb_atome1 != nb_atome4:
        print("\033[31mErreur: Les PDB ne comportent %s\
pas le meme nombre d'atomes apres selection\033[0m")
        sys.exit()
    vegan = importr("vegan")
    matricedif1 = []
    matricedif2 = []
    for i in range(nb_atome1*nb_atome1):
        matricedif1.append(matrice1[i] - matrice2[i]) 
        matricedif2.append(matrice3[i] - matrice4[i])
    v1 = robjects.FloatVector(matricedif1)
    v2 = robjects.FloatVector(matricedif2)
    m1 = robjects.r['matrix'](v1, nrow = nb_atome1)
    m2 = robjects.r['matrix'](v2, nrow = nb_atome1)
    coeff = vegan.mantel(m1, m2)
    return coeff[2][0]


def CoeffCorrelationR(Omatrice1, Omatrice2, Omatrice3, Omatrice4):
    '''Calcul le coeff de correlation entre deux mouvements (-> 4 matrices)'''
    try:
        from rpy2.robjects.packages import importr
        import rpy2.robjects as robjects
    except ImportError:
        print("\033[31mFonction : CoeffCorrelationCompletR")
        print("Module rpy2 requis !\033[0m")
        sys.exit()
    (matrice1, nb_atome1) = Omatrice1
    (matrice2, nb_atome2) = Omatrice2
    (matrice3, nb_atome3) = Omatrice3
    (matrice4, nb_atome4) = Omatrice4
    if nb_atome1 != nb_atome2 or nb_atome1 != nb_atome3 or \
       nb_atome1 != nb_atome4:
        print("\033[31mErreur: Les PDB ne comportent %s pas le meme nombre \
d'atomes apres selection\033[0m")
        sys.exit()
    vegan = importr("vegan")
    matricedif1 = []
    matricedif2 = []
    for i in range(nb_atome1*nb_atome1):
        matricedif1.append(matrice1[i] - matrice2[i]) 
        matricedif2.append(matrice3[i] - matrice4[i])
    v1 = robjects.FloatVector(matricedif1)
    v2 = robjects.FloatVector(matricedif2)
    m1 = robjects.r['matrix'](v1, nrow = nb_atome1)
    m2 = robjects.r['matrix'](v2, nrow = nb_atome1)
    coeff = vegan.mantel(m1, m2)
    return coeff[2][0]


def CoeffCorrelation(matriceDep1, matriceDep2):
    '''Calcul le coeff de correlation entre deux mouvements'''
    res1 = 0
    res2 = 0
    res3 = 0
    longueur = len(matriceDep1)
    moyenne1 = sum(matriceDep1) / longueur
    moyenne2 = sum(matriceDep2) / longueur
    for i in range(longueur):
        res1 += ((matriceDep1[i] - moyenne1) * (matriceDep2[i] - moyenne2))
        res2 += pow(matriceDep1[i] - moyenne1, 2)
        res3 += pow(matriceDep2[i] - moyenne2, 2)
    res_f = res1 / ( sqrt(res2) * sqrt(res3) )
    return res_f


##############################################################################
###### MOTION COMBINATION   ##################################################
##############################################################################
def MotionCombination(PDB1, PDB2, PDB3, PDB4, PDB5, atom_select):
    selection = PDB1.CheckSel(atom_select, "a")
    liste_atome1 = []
    for a in PDB1.liste_atome:
        test = eval(selection)
        if test:
            liste_atome1.append(a)
    selection = PDB2.CheckSel(atom_select, "a")
    liste_atome2 = []
    for a in PDB2.liste_atome:
        test = eval(selection)
        if test:
            liste_atome2.append(a)
    selection = PDB1.CheckSel(atom_select, "a")
    liste_atome3 = []
    for a in PDB3.liste_atome:
        test = eval(selection)
        if test:
            liste_atome3.append(a)
    selection = PDB1.CheckSel(atom_select, "a")
    liste_atome4 = []
    for a in PDB4.liste_atome:
        test = eval(selection)
        if test:
            liste_atome4.append(a)
    for i in range(len(liste_atome1)):
        xa1 = liste_atome1[i].x
        xb1 = liste_atome2[i].x
        xa2 = liste_atome3[i].x
        xb2 = liste_atome4[i].x
        ya1 = liste_atome1[i].y
        yb1 = liste_atome2[i].y
        ya2 = liste_atome3[i].y
        yb2 = liste_atome4[i].y
        za1 = liste_atome1[i].z
        zb1 = liste_atome2[i].z
        za2 = liste_atome3[i].z
        zb2 = liste_atome4[i].z
        componentx = (xa1 - xb1) + (xa2 - xb2)
        componenty = (ya1 - yb1) + (ya2 - yb2)
        componentz = (za1 - zb1) + (za2 - zb2)
        PDB5.liste_atome[i].x += componentx/2
        PDB5.liste_atome[i].y += componenty/2
        PDB5.liste_atome[i].z += componentz/2


##############################################################################
###### BFACTOR REMPLACEMENT ##################################################
##############################################################################
def Bfactor(PDB, bfactorFile, atom_select):
    '''Assigne les Bfactor contenus dans un fichier sur un PDB'''
    try:
        fichier_in_b = open(bfactorFile, 'r')
    except IOError:
        print("\033[31mErreur: Impossible d'ouvrir le fichier %s\033[0m"\
              %(bfactorFile))
        sys.exit()
    else:
        b_factor = []
        for ligne in fichier_in_b.readlines():
            b_factor.append(float(ligne.split()[-1]))
        fichier_in_b.close()
    selection = PDB.CheckSel(atom_select, "a")
    compteur_atm = 0
    for a in PDB.liste_atome:
        test = eval(selection)
        if test:
            a.bfactor = b_factor[compteur_atm]
            compteur_atm += 1


##############################################################################
###### ECHANGE NOM ATOME   ###################################################
##############################################################################
def EchangeNomAtome(PDB1, PDB2):
    '''Met le nom des atomes du PDB1 dans le PDB2'''
    if len(PDB1.liste_atome) != len(PDB2.liste_atome):
        print("Erreur: Les deux PDB doivent comporter le meme nombre d'atome")
    index = 0
    for a1 in PDB1.liste_atome:
        PDB2.liste_atome[index].type_atome = a1.type_atome
        index += 1
    return PDB2


##############################################################################
###### CALUL L'ANGLE DE DEUX VECTEURS ENTRE DEUX PDB #########################
##############################################################################
def AngleHeliceRef(PDBRef, PDB2, atom_select):
    '''Calcul l'angle entre deux vecteurs (2 atomes) de deux PDB differents'''
    testpdbref = isinstance(PDBRef, pdb)
    testpdb = isinstance(PDB2, pdb)
    testsel = isinstance(atom_select, str)
    if not testpdbref or not testpdb or not testsel:
        print("\033[31mFonction : AngleHeliceRef")
        print("Format argument(s) incorrect(s) !\033[0m")
        sys.exit()
    selection1 = PDBRef.CheckSel(atom_select, "a1")
    selection2 = PDB2.CheckSel(atom_select, "a2")
    liste1 = []
    liste2 = []
    for a1 in PDBRef.liste_atome:
        test1 = eval(selection1)
        if test1:
            liste1.append(a1)
    for a2 in PDB2.liste_atome:
        test2 = eval(selection2)
        if test2:
            liste2.append(a2)
    xa_tar = liste2[0].x
    ya_tar = liste2[0].y
    za_tar = liste2[0].z
    xb_tar = liste2[-1].x
    yb_tar = liste2[-1].y
    zb_tar = liste2[-1].z
    xa_ref = liste1[0].x
    ya_ref = liste1[0].y
    za_ref = liste1[0].z
    xb_ref = liste1[-1].x
    yb_ref = liste1[-1].y
    zb_ref = liste1[-1].z
    xa = xb_tar - xa_tar
    ya = yb_tar - ya_tar
    za = zb_tar - za_tar
    xb = xb_ref - xa_ref
    yb = yb_ref - ya_ref
    zb = zb_ref - za_ref
    xb_rot = xb_tar - xa_ref
    yb_rot = yb_tar - ya_ref
    zb_rot = zb_tar - za_ref
    if xa + ya + zb >= xb + yb + zb:
        angle = -180 * (acos((xa*xb + ya*yb + za*zb) / \
               (sqrt((xa*xa + ya*ya + za*za) * \
               (xb*xb + yb*yb + zb*zb))))) / pi
        angle_rot = -180 * (acos((xa*xb_rot + ya*yb_rot + za*zb_rot) / \
                    (sqrt((xa*xa + ya*ya + za*za) * \
                    (xb_rot*xb_rot + yb_rot*yb_rot + zb_rot*zb_rot))))) / \
                    pi
    elif xa + ya + zb < xb + yb + zb:
        angle = 180 * (acos((xa*xb + ya*yb + za*zb) / \
                (sqrt((xa*xa + ya*ya + za*za) * (xb*xb + yb*yb + zb*zb))))) / \
                pi
        angle_rot = 180 * (acos((xa*xb_rot + ya*yb_rot + za*zb_rot) / \
                    (sqrt((xa*xa + ya*ya + za*za) * \
                    (xb_rot*xb_rot + yb_rot*yb_rot + zb_rot*zb_rot))))) / pi
    return [angle, angle_rot]


def SelfCrossCorrelationMatrixZeros(premier_mode, dernier_mode, atom_select,
                                    path_origin, matrice_out):
    mode = premier_mode
    matrice_liste = []
    while mode <= dernier_mode:
        matriceres = []
        PDB1 = pdb(path_origin + "/mini-mode-" + str(mode) + \
                   "/conformation_1.pdb", force=True)
        PDB2 = pdb(path_origin + "/mini-mode-" + str(mode) + \
                   "/conformation_-1.pdb", force=True)
        matrice1 = PDB1.DistanceMatrixR(atom_select)
        matrice2 = PDB2.DistanceMatrixR(atom_select)
        for i in range(len(matrice1[0])):
            matriceres.append(matrice1[0][i] - matrice2[0][i])
        matrice_liste.append(matriceres)
        mode += 1
    fo = open(matrice_out, 'w')
    for i in range(dernier_mode - premier_mode + 1):
        fo.write("   0.000000" * i)
        for j in range(dernier_mode - premier_mode - i + 1):
            fo.write(("%1.6f" %(CoeffCorrelation(matrice_liste[i], matrice_liste[j + i]))).rjust(11))
            print("Coefficients Modes : %d - %d effectue" %(i + premier_mode,
                                                            j + i + \
                                                            premier_mode))
        fo.write("\n")
    fo.close()


def SelfCrossCorrelationMatrix(premier_mode, dernier_mode, atom_select,
                               path_origin, matrice_out):
    mode = premier_mode
    matrice_liste = []
    while mode <= dernier_mode:
        matriceres = []
        PDB1 = pdb(path_origin + "/mini-mode-" + str(mode) + \
                   "/conformation_1.pdb", force=True)
        PDB2 = pdb(path_origin + "/mini-mode-" + str(mode) + \
                   "/conformation_-1.pdb", force=True)
        matrice1 = PDB1.DistanceMatrixR(atom_select)
        matrice2 = PDB2.DistanceMatrixR(atom_select)
        for i in range(len(matrice1[0])):
            matriceres.append(matrice1[0][i] - matrice2[0][i])
        matrice_liste.append(matriceres)
        mode += 1
    mat = []
    for i in range(dernier_mode - premier_mode + 1):
        l = []
        for j in range(dernier_mode - premier_mode + 1):
            l.append(0.0)
        mat.append(l)
    for i in range(dernier_mode - premier_mode + 1):
        for j in range(dernier_mode - premier_mode - i + 1):
            coeff = CoeffCorrelation(matrice_liste[i], matrice_liste[j + i])
            mat[i][j + i] = coeff
            mat[j + i ][i] = coeff
            print("Coefficients Modes : %d - %d effectue" %(i + premier_mode,
                                                            j + i + \
                                                            premier_mode))
    fo = open(matrice_out, 'w')
    for ligne in mat:
        for valeur in ligne:
            fo.write(("%1.6f" % valeur).rjust(11))
        fo.write("\n")
    fo.close()


def CrossCorrelationMatrix(premier_mode1, dernier_mode1, premier_mode2,
                           dernier_mode2, atom_select1, path_origin1,
                           atom_select2, path_origin2, matrice_out):
    mode1 = premier_mode1
    mode2 = premier_mode2
    matrice_liste1 = []
    matrice_liste2 = []
    while mode1 <= dernier_mode1:
        matriceres1 = []
        PDB1 = pdb(path_origin1 + "/mini-mode-" + str(mode1) + \
                   "/conformation_1.pdb", force=True)
        PDB2 = pdb(path_origin1 + "/mini-mode-" + str(mode1) + \
                   "/conformation_-1.pdb", force=True)
        matrice1 = PDB1.DistanceMatrixR(atom_select1)
        matrice2 = PDB2.DistanceMatrixR(atom_select1)
        del PDB1, PDB2
        if matrice1[1] != matrice2[1]:
            print("Erreur: Les PDB ne comportent pas le meme nombre d'atome \
apres selection ! (", matrice1[1], matrice2[1], ")")
            sys.exit()
        taille = len(matrice1[0])
        nb_atom = matrice1[1]
        for i in range(taille):
            matriceres1.append(matrice1[0][i] - matrice2[0][i])
        del matrice1, matrice2
        matrice_liste1.append(matriceres1)
        mode1 += 1
    while mode2 <= dernier_mode2:
        matriceres2 = []
        PDB3 = pdb(path_origin2 + "/mini-mode-" + str(mode2) + \
                   "/conformation_1.pdb", force=True)
        PDB4 = pdb(path_origin2 + "/mini-mode-" + str(mode2) + \
                   "/conformation_-1.pdb", force=True)
        matrice3 = PDB3.DistanceMatrixR(atom_select2)
        matrice4 = PDB4.DistanceMatrixR(atom_select2)
        del PDB3, PDB4
        if nb_atom != matrice3[1] or nb_atom != matrice4[1]:
            print("Erreur: Les PDB ne comportent pas le meme nombre d'atome \
apres selection ! (", nb_atom, matrice3[1], matrice4[1], ")")
            sys.exit()
        for i in range(taille):
            matriceres2.append(matrice3[0][i] - matrice4[0][i])
        del matrice3, matrice4
        matrice_liste2.append(matriceres2)
        mode2 += 1
    fo = open(matrice_out, 'w')
    for i in range(dernier_mode1 - premier_mode1 + 1):
        for j in range(dernier_mode2 - premier_mode2 + 1):
            fo.write(("%1.6f" %(CoeffCorrelation(matrice_liste1[i], matrice_liste2[j]))).rjust(11))
            print("Coefficients Modes : %d - %d effectue" %(i + premier_mode1,
                                                            j + premier_mode2))
        fo.write("\n")
    fo.close()
    del matrice_liste1, matrice_liste2


def Distance(PDB1, atom_select1, PDB2, atom_select2):
    '''Calcul la distance entre deux selections d'atomes (centre de masse)'''
    testpdb1 = isinstance(PDB1, pdb)
    testpdb2 = isinstance(PDB2, pdb)
    testsel1 = isinstance(atom_select1, str)
    testsel2 = isinstance(atom_select2, str)
    if not testpdb1 or not testpdb2 or not testsel1 or not testsel2:
        print("\033[31mFonction : Distance")
        print("Format argument(s) incorrect(s) !\033[0m")
        sys.exit()
    cm1 = PDB1.CentreDeMasse(atom_select1)
    cm2 = PDB2.CentreDeMasse(atom_select2)
    distance = sqrt(pow(cm2[0] - cm1[0], 2) + pow(cm2[1] - cm1[1], 2) + \
                    pow(cm2[2] - cm1[2], 2))
    return distance


def TMD(PDBi, atom_select1, PDBt, atom_select2):
    '''Prepare un fichier de TMD avec un PDB ini et target et une selection'''
    # Remplace les coordonnees des atomes selectionnes par atom_select1 de PDBi
    # par ceux selectionnes par atom_select2 de PDBt
    PDBres = PDBi.Copy("all")
    PDBres.ModifOccupancy(0.0, "all")
    Lindex1 = []
    Lindex2 = []
    selection1 = PDBi.CheckSel(atom_select1, "a")
    selection2 = PDBt.CheckSel(atom_select2, "a")
    i = 0
    for a in PDBi.liste_atome:
        if eval(selection1):
            Lindex1.append(i)
        i += 1
    i = 0
    for a in PDBt.liste_atome:
        if eval(selection2):
            Lindex2.append(i)
        i += 1
    taille = len(Lindex1)
    if taille != len(Lindex2):
        print("Erreur: les selections n'ont pas le meme nombre d'atome !")
        sys.exit()
    for i in range(taille):
        PDBres.liste_atome[Lindex1[i]].x = PDBt.liste_atome[Lindex2[i]].x
        PDBres.liste_atome[Lindex1[i]].y = PDBt.liste_atome[Lindex2[i]].y
        PDBres.liste_atome[Lindex1[i]].z = PDBt.liste_atome[Lindex2[i]].z
        PDBres.liste_atome[Lindex1[i]].occupancy = 1.0
    return PDBres


def pdb2vec(PDBref, PDBini, PDBfinal, normalization=False):
    '''Genere un objet pdb() copie de PDBref avec comme coordonnees vecteur
    resultant du deplacement entre PDBini et PDBfinal, met le Bfactor a 1 pour
    les atomes selectionnes et 0 pour les autres. Les atomes selectionnes sont
    ceux communs entre le mouvement (PDBfinal - PDBini) et le PDBref, le script
    regarde type_res, num_res, chain et type_atome. Le vecteur est normalise.
    '''
    # Remplace les coordonnees des atomes selectionnes par atom_select1 de PDBi
    # par ceux selectionnes par atom_select2 de PDBt
    PDBres = PDBref.Copy("all")
    PDBvec = PDBini.Copy("all")
    vec = []
    minx = PDBini.liste_atome[0].x - PDBfinal.liste_atome[0].x
    miny = PDBini.liste_atome[0].y - PDBfinal.liste_atome[0].y
    minz = PDBini.liste_atome[0].z - PDBfinal.liste_atome[0].z
    maxx = PDBini.liste_atome[0].x - PDBfinal.liste_atome[0].x
    maxy = PDBini.liste_atome[0].y - PDBfinal.liste_atome[0].y
    maxz = PDBini.liste_atome[0].z - PDBfinal.liste_atome[0].z
    for i in range(len(PDBini.liste_atome)):
        a1 = PDBini.liste_atome[i]
        a2 = PDBfinal.liste_atome[i]
        x, y, z = a2.x - a1.x, a2.y - a1.y, a2.z - a1.z
        if normalization:
            norm = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
            vec.append([x / norm, y / norm, z / norm])
        else:
            vec.append([x, y, z])
    for i in range(len(PDBvec.liste_atome)):
        PDBvec.liste_atome[i].x = vec[i][0]
        PDBvec.liste_atome[i].y = vec[i][1]
        PDBvec.liste_atome[i].z = vec[i][2]
    for a1 in PDBres.liste_atome:
        found = False
        for a2 in PDBvec.liste_atome:
            if a1.type_res == a2.type_res and a1.num_res == a2.num_res and \
               a1.chain == a2.chain and a1.type_atome == a2.type_atome:
                found = True
                a1.x, a1.y, a1.z = a2.x, a2.y, a2.z
                break
        if not found:
            a1.bfactor = 0.00
            a1.x, a1.y, a1.z = 0.0, 0.0, 0.0
        else:
            a1.bfactor = 1.00
        a1.occupancy = 0.00
    return PDBres


def AjoutOxygenHem(PDB, PDBref):
    for a in PDBref.liste_atome:
        if a.type_res == "HEO" and a.type_atome == "FE":
            FeRef = a
        elif a.type_res == "HEO" and a.type_atome == "OE":
            ORef = a
        elif a.type_res == "HEO" and a.type_atome == "NA":
            NaRef = a
        elif a.type_res == "HEO" and a.type_atome == "NB":
            NbRef = a
    for a in PDB.liste_atome:
        if a.type_res == "HEM" and a.type_atome == "FE":
            Fe = a
        if a.type_res == "HEM" and a.type_atome == "NA":
            Na = a
        if a.type_res == "HEM" and a.type_atome == "NB":
            Nb = a
    xO = (ORef.x - FeRef.x) + Fe.x
    yO = (ORef.y - FeRef.y) + Fe.y
    zO = (ORef.z - FeRef.z) + Fe.z
    NOrefx = ORef.x - NaRef.x
    NOrefy = ORef.y - NaRef.y
    NOrefz = ORef.z - NaRef.z
    NOx = xO - Na.x
    NOy = yO - Na.y
    NOz = zO - Na.z
    xO = NOx - NOrefx
    yO = NOy - NOrefy
    zO = NOz - NOrefz
    newO = atome("ATOM", 1, "OE", "HEO", "A", 501, xO, yO, zO, 1.00, 0.00, 
                 "HEO ", "   ", 1.00)
    PDB.AjoutAtm(newO)
    PDB.RenumAtm(1, "all")
    return PDB


def AjoutOxygenHem2(PDB):
    for a in PDB.liste_atome:
        if a.type_res == "HEME" and a.type_atome == "FE":
            Fe = a
        if a.num_res == 435 and a.type_atome == "SG":
            Sg = a
    xO = 0.6867 * (Fe.x - Sg.x) + Fe.x
    yO = 0.6867 * (Fe.y - Sg.y) + Fe.y
    zO = 0.6867 * (Fe.z - Sg.z) + Fe.z
    newO = atome("ATOM", 1, "OE", "HEO", "A", 501, xO, yO, zO, 1.00, 0.00, 
                 "HEO ", "   ", 1.00)
    PDB.AjoutAtm(newO)
    PDB.RenumAtm(1, "all")
    return PDB


def TransformPDB(PDB, vector, atom_select):
    #vector is a list of length N (number of selected atoms) containing
    #list of three values (x, y and z)
    newPDB = PDB.Copy("all")
    selection = newPDB.CheckSel(atom_select, "a")
    selection = selection + " and a.desc == \"ATOM\""
    i=0
    for a in newPDB.liste_atome:
        test = eval(selection)
        if test:
            a.x = a.x - vector[i][0]
            a.y = a.y - vector[i][1]
            a.z = a.z - vector[i][2]
            i += 1
    return newPDB


def COMFromListAtm(listAtom):
    x, y, z, masse = 0.0, 0.0, 0.0, 0.0
    for a in listAtom:
        x += a.x * a.mass
        y += a.y * a.mass
        z += a.z * a.mass
        masse += a.mass
    return x / masse, y / masse, z / masse


def CenterFromListAtm(listAtom):
    x, y, z = 0.0, 0.0, 0.0
    for a in listAtom:
        x += a.x
        y += a.y
        z += a.z
    return x / len(listAtom), y / len(listAtom), z / len(listAtom)


def AllAtom2CG(PDB, atom_select):
    cgPDB = pdb()
    selection = PDB.CheckSel(atom_select + " and protein", "a") + " and a.desc == \"ATOM\""
    num_atome = 1
    previousAtm = atome("ATOM", 1, "X", "X", "X", -10, 0.0, 0.0, 0.0, 0.0, 0.0,
                        "X", "X")
    listAtomBB, listAtomSC = [], []
    for a in PDB.liste_atome:
        test = eval(selection)
        if test:
            if previousAtm.num_res != a.num_res or previousAtm.type_res != a.type_res:
                if len(listAtomBB) > 1:
                    xbb, ybb, zbb = CenterFromListAtm(listAtomBB)
                    bb = atome("ATOM", num_atome, "BB", listAtomBB[0].type_res,
                           listAtomBB[0].chain, listAtomBB[0].num_res, xbb, ybb,
                           zbb, 1.0, 2.0, listAtomBB[0].segname, "X")
                    cgPDB.AjoutAtm(bb)
                    num_atome += 1
                elif len(listAtomBB) == 1:
                    bb = atome("ATOM", num_atome, "BB", listAtomBB[0].type_res,
                               listAtomBB[0].chain, listAtomBB[0].num_res,
                               listAtomBB[0].x, listAtomBB[0].y,
                               listAtomBB[0].z, 1.0, 2.0, listAtomBB[0].segname,
                               "X")
                    cgPDB.AjoutAtm(bb)
                    num_atome += 1
                if len(listAtomSC) > 1:
                    xsc, ysc, zsc = CenterFromListAtm(listAtomSC)
                    sc = atome("ATOM", num_atome, "SC", listAtomSC[0].type_res,
                           listAtomSC[0].chain, listAtomSC[0].num_res, xsc, ysc,
                           zsc, 1.0, a.cgRadii[previousAtm.type_res],
                           listAtomSC[0].segname, "X")
                    cgPDB.AjoutAtm(sc)
                    num_atome += 1
                elif len(listAtomSC) == 1:
                    sc = atome("ATOM", num_atome, "SC", listAtomSC[0].type_res,
                               listAtomSC[0].chain, listAtomSC[0].num_res,
                               listAtomSC[0].x, listAtomSC[0].y,
                               listAtomSC[0].z, 1.0,
                               a.cgRadii[previousAtm.type_res],
                               listAtomSC[0].segname, "X")
                    cgPDB.AjoutAtm(sc)
                    num_atome += 1
                listAtomBB, listAtomSC = [], []
            if a.type_atome in ["CA", "N", "C", "O", "HN"]:
                listAtomBB.append(atome("ATOM", num_atome, "BB", a.type_res,
                                        a.chain, a.num_res, a.x, a.y, a.z,
                                        1.0, 1.0, a.segname, a.element))
            elif a.type_res != "GLY":
                listAtomSC.append(atome("ATOM", num_atome, "SC", a.type_res,
                                        a.chain, a.num_res, a.x, a.y, a.z,
                                        1.0, 1.0, a.segname, a.element))
            previousAtm = a
    connect = ""
    previousAtm = atome("ATOM", 1, "X", "X", "X", -10, 0.0, 0.0, 0.0, 0.0, 0.0,
                        "X", "X")
    for i in range(len(cgPDB.liste_atome)):
        try:
            a = cgPDB.liste_atome[i]
            amoins1 = cgPDB.liste_atome[i-1]
            amoins2 = cgPDB.liste_atome[i-2]
        except IndexError:
            pass
        else:
            if amoins1.type_atome == "BB" and a.type_atome == "SC" and a.num_res == amoins1.num_res:
                connect += "CONECT%5d%5d\n" %(amoins1.num_atome, a.num_atome)
            if amoins1.type_atome == "BB" and a.type_atome == "BB" and (a.num_res-1) == amoins1.num_res:
                connect += "CONECT%5d%5d\n" %(amoins1.num_atome, a.num_atome)
            if amoins2.type_atome == "BB" and a.type_atome == "BB" and (a.num_res-1) == amoins2.num_res:
                connect += "CONECT%5d%5d\n" %(amoins2.num_atome, a.num_atome)
        previousAtm = cgPDB.liste_atome[i]
    cgPDB.connect = connect
    return cgPDB


def addVec2PDB(PDB, vecFile):
    newPDB = PDB.Copy("all")
    fi = open(vecFile, 'r')
    for line in fi.readlines():
        items = line.split()
        atomId = int(items[1])
        selection = newPDB.CheckSel("index %d" %(atomId), "a")
        for a in newPDB.liste_atome:
            test = eval(selection)
            if test:
                a.x = a.x + float(items[2])
                a.y = a.y + float(items[3])
                a.z = a.z + float(items[4])
                break
    fi.close()
    return newPDB


def subtractVec2PDB(PDB, vecFile):
    newPDB = PDB.Copy("all")
    fi = open(vecFile, 'r')
    for line in fi.readlines():
        items = line.split()
        atomId = int(items[1])
        selection = newPDB.CheckSel("index %d" %(atomId), "a")
        for a in newPDB.liste_atome:
            test = eval(selection)
            if test:
                a.x = a.x - float(items[2])
                a.y = a.y - float(items[3])
                a.z = a.z - float(items[4])
                break
    fi.close()
    return newPDB

def distance(a1, a2):
    dist = sqrt(pow(a2.x - a1.x, 2) + pow(a2.y - a1.y, 2) + pow(a2.z - a1.z, 2))
    return dist


if __name__ == "__main__":
    ############################################################################
    ###### DEBUT PARTIE SCRIPT #################################################
    ############################################################################
    PDB = pdb("/home/ali/RELAXIN/alphafold/start_cg.pdb",verbose=True, force=True)
    PDB.RenumRes(20, "chain A")
    PDB.ExportPDB("/home/ali/RELAXIN/alphafold/start_cg_renum.pdb", "all", TER=True)

    ############################################################################
    ###### FIN PARTIE SCRIPT ###################################################
    ############################################################################
    #ressources = resource.getrusage(resource.RUSAGE_SELF)
    #print("=================================")
    #print("Temps utilisateur:    %.3f sec" %ressources[0])
    #print("Temps systeme:        %.3f sec" %ressources[1])
    #print("Memoire Vive Maximum: %.2f Mo" %(ressources[2] / 1024))
    #print("=================================")
