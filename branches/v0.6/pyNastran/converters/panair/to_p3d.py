from panairGrid import PanairGrid

def run():
    infileName = 'HSCT.inp'
    model = PanairGrid(infileName, log=None, debug=True)
    model.read_panair()

    p3d_name = 'HSCT.plt'
    model.write_plot3d(p3d_name)

if __name__ == '__main__':
    run()
