from sdpp import read

def convert(file):
  print("Reading File")
  data,meta = read(file)

  print("Writing ROA File")
  with open(file.rstrip(".dat")+".roa",'w') as f:
    f.write("\nTYPE ROA\n")
    f.write("UF doublesidedstrip adc\n")
    for event_id in range(len(data)):
      det = 0    # arbitrary (only one detector)
      strip = 1  # arbitrary (only one strip)
      side = 'n' # arbitrary (only one side)
 
      f.write("\nSE\n")
      f.write("ID {}\n".format(event_id))
      f.write("CL {:f}\n".format(data['time'][event_id]))
      f.write("UH {} {} {} {}\n".format(det, strip, side, int(data['amplitude'][event_id])))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Convert discrete preamp raw data files to Melinator ROA files')
    parser.add_argument('file', default=None, nargs='?')
    args = parser.parse_args()

    convert(args.file)
