
#         Makefile for SOSIE 
#         ******************

# Version: trunk


# Create a symbolic link from macro/make.macro_* of your choice to ./make.macro !!!

include make.macro

NCDF_LIB=$(NETCDF_DIR)/lib
NCDF_INC=$(NETCDF_DIR)/include


#========================================================
#          You should not change anything below
#========================================================

LIB_SOSIE = lib/libsosie.a

# LIBRARY :
LIB_CDF = -L$(NCDF_LIB) $(L_NCDF)
LIB     = -L./lib -lsosie $(LIB_CDF)


OBJ = obj/io_ezcdf.o \
      obj/mod_conf.o \
      obj/mod_init.o \
      obj/mod_manip.o \
      obj/mod_grids.o \
      obj/mod_drown.o \
      obj/mod_akima_2d.o \
      obj/mod_bilin_2d.o \
      obj/mod_akima_1d.o \
      obj/mod_scoord.o \
      obj/mod_interp.o \
      obj/mod_nemotools.o \
      obj/mod_strftime.o \
      obj/mod_poly.o


OBJ_I2GT = obj/io_ezcdf.o \
           obj/mod_conf.o \
           obj/mod_manip.o \
           obj/mod_drown.o \
           obj/mod_bilin_2d.o \
           obj/mod_poly.o


# Modules to install in $INSTALL_DIR/include :
MOD_INST= mod/io_ezcdf.mod \
          mod/mod_akima_2d.mod \
          mod/mod_drown.mod

all: bin/sosie.x bin/corr_vect.x bin/mask_drown_field.x bin/interp_to_ground_track.x bin/ij_from_lon_lat.x

gt: bin/interp_to_ground_track.x

test: bin/test_stuffs.x bin/test_poly.x

bin/sosie.x: src/sosie.f90 $(LIB_SOSIE)
	@mkdir -p bin
	$(FC) $(FF) src/sosie.f90 -o bin/sosie.x $(LIB)

bin/corr_vect.x: src/corr_vect.f90 $(LIB_SOSIE)
	$(FC) $(FF) src/corr_vect.f90 -o bin/corr_vect.x $(LIB)

bin/test_stuffs.x: src/test_stuffs.f90 $(LIB_SOSIE)
	$(FC) $(FF) src/test_stuffs.f90 -o bin/test_stuffs.x $(LIB)

bin/test_poly.x: src/test_poly.f90 $(LIB_SOSIE)
	$(FC) $(FF) src/test_poly.f90 -o bin/test_poly.x $(LIB)



bin/interp_to_ground_track.x: src/interp_to_ground_track.f90 $(OBJ_I2GT)
	@mkdir -p bin
	$(FC) $(FF) $(OBJ_I2GT) src/interp_to_ground_track.f90 -o bin/interp_to_ground_track.x $(LIB_CDF)

bin/ij_from_lon_lat.x: src/ij_from_lon_lat.f90 obj/io_ezcdf.o obj/mod_manip.o
	@mkdir -p bin
	$(FC) $(FF) obj/io_ezcdf.o obj/mod_manip.o src/ij_from_lon_lat.f90 -o bin/ij_from_lon_lat.x $(LIB_CDF)




$(LIB_SOSIE): $(OBJ)
	@echo
	@mkdir -p lib bin
	ar -rv $(LIB_SOSIE) $(OBJ)
	ranlib $(LIB_SOSIE)
	@echo



# Objects :
# ---------

obj/io_ezcdf.o: src/io_ezcdf.f90
	@mkdir -p obj
	@mkdir -p mod
	$(FC) $(FF) -I$(NCDF_INC) -c src/io_ezcdf.f90 -o obj/io_ezcdf.o

obj/mod_conf.o: src/mod_conf.f90 obj/io_ezcdf.o
	@mkdir -p obj
	@mkdir -p mod
	$(FC) $(FF) -c src/mod_conf.f90 -o obj/mod_conf.o

obj/mod_init.o: src/mod_init.f90 obj/mod_conf.o obj/mod_scoord.o
	$(FC) $(FF) -c src/mod_init.f90 -o obj/mod_init.o



obj/mod_grids.o: src/mod_grids.f90 obj/mod_conf.o obj/io_ezcdf.o obj/mod_manip.o
	$(FC) $(FF) -c src/mod_grids.f90 -o obj/mod_grids.o

obj/mod_interp.o: src/mod_interp.f90 obj/mod_nemotools.o
	$(FC) $(FF) -c src/mod_interp.f90 -o obj/mod_interp.o

obj/mod_manip.o: src/mod_manip.f90
	$(FC) $(FF) -c src/mod_manip.f90 -o obj/mod_manip.o

obj/mod_drown.o: src/mod_drown.f90
	$(FC) $(FF) -c src/mod_drown.f90 -o obj/mod_drown.o

obj/mod_akima_2d.o: src/mod_akima_2d.f90
	$(FC) $(FF) -c src/mod_akima_2d.f90 -o obj/mod_akima_2d.o

obj/mod_bilin_2d.o: src/mod_bilin_2d.f90 obj/io_ezcdf.o obj/mod_conf.o obj/mod_manip.o obj/mod_poly.o
	$(FC) $(FF) -c src/mod_bilin_2d.f90 -o obj/mod_bilin_2d.o

obj/mod_akima_1d.o: src/mod_akima_1d.f90
	$(FC) $(FF) -c src/mod_akima_1d.f90 -o obj/mod_akima_1d.o

obj/mod_scoord.o: src/mod_scoord.f90
	$(FC) $(FF) -c src/mod_scoord.f90 -o obj/mod_scoord.o

obj/mod_nemotools.o: src/mod_nemotools.f90
	$(FC) $(FF) -c src/mod_nemotools.f90 -o obj/mod_nemotools.o

obj/mod_strftime.o: src/mod_strftime.f90
	$(FC) $(FF) -c src/mod_strftime.f90 -o obj/mod_strftime.o

obj/mod_poly.o: src/mod_poly.f90 obj/io_ezcdf.o
	$(FC) $(FF) -c src/mod_poly.f90 -o obj/mod_poly.o




bin/mask_drown_field.x: src/mask_drown_field.f90 $(LIB_SOSIE)
	$(FC) $(FF) -o bin/mask_drown_field.x src/mask_drown_field.f90 $(LIB)


install: all
	@mkdir -p $(INSTALL_DIR)/bin $(INSTALL_DIR)/lib $(INSTALL_DIR)/include
	cp bin/*.x      $(INSTALL_DIR)/bin/
	cp $(LIB_SOSIE) $(INSTALL_DIR)/lib/
	cp $(MOD_INST)  $(INSTALL_DIR)/include/


uninstall:
	rm -f $(INSTALL_DIR)/bin/*.x
	rm -f $(INSTALL_DIR)/lib/$(LIB_SOSIE) $(INSTALL_DIR)/bin/apply_zon_corr.x
	rm -f $(INSTALL_DIR)/include/io_ezcdf.mod


clean:
	rm -f bin/* $(LIB_SOSIE) *.out mod/*.mod *.x *~ *\# src/*~ src/*\# *.log
	rm -rf bin mod lib obj
	rm -f examples/*.nc4





