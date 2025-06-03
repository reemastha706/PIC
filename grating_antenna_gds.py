import gdsfactory as gf
import random

def waveguide(l, w, ly):
  xs = gf.cross_section.strip(width = w, layer = [ly, 0])
  return gf.components.straight(length = l, cross_section=xs)

def waveguide_antenna(wg, h):
  
  w0 = 1  #substrate
  w1 = 3  # sio2
  w2 = wg  #si waeguide  500 nm = 0.5
  n_periods, period, duty_cycle, grating_height = 5, 0.63, 0.5 ,h
  idx = random.random()
  c = gf.Component()

  l = period * n_periods+2
  
  l0 = c << waveguide(l, w0, 1)
  l0.move((0, w0/2))
  l1 = c << waveguide(l, w1, 2)
  l1.move((0,w0+w1/2))
  l2 = c << waveguide(l, w2, 3)
  l2.move((0, w0+ w1 + w2/2))

  for i in range(n_periods):
        tooth = c << waveguide(
        l = period * duty_cycle,
        w=grating_height,
        ly=4
        )
        x_position = c.xsize/2 -l/2 + 1 + i * period + (period * (1 - duty_cycle) / 2)
        tooth.move(( x_position, w0+w1+w2+grating_height/2))
        
  c.flatten()


  gdsII_file= "path"+str(idx)+".gds"

  c.write_gds(gdsII_file)
  return c.xsize, c.ysize, gdsII_file