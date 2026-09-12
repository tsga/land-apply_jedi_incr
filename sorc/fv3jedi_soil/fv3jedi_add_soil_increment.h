#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

#include "oops/base/Variables.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "soil_increments_cpp_interface.h"

namespace landincr {
  /**
   * AddLandIncrement Class Implementation
   *
   */

  class AddLandIncrement : public oops::Application {
   public:
    explicit AddLandIncrement(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}

    virtual ~AddLandIncrement() {}

    static const std::string classname() {return "landincr::AddLandIncrement";}

    int execute(const eckit::Configuration & fullConfig) const override {

      int myrank = this->getComm().rank();

      // We assume both state and increment are in the same geometry
      // const ijedi::Geometry<fv3jedi::Traits> geom_(eckit::LocalConfiguration(fullConfig, "geometry"),this->getComm());
      fv3jedi::Geometry geom_(eckit::LocalConfiguration(fullConfig, "geometry"),this->getComm());
      // oops::Log::info() << "geom ny "<<geom_.npy() << " nx " << geom_.npx() << std::endl;

      // Read state
      // ijedi::State<ijedi::Traits> xx(geom_, eckit::LocalConfiguration(fullConfig, "background state"));
      fv3jedi::State xx(geom_, eckit::LocalConfiguration(fullConfig, "background state"));
      oops::Log::test() << "Background state: " << xx << std::endl;

      atlas::FieldSet bkg_fs;
      xx.toFieldSet(bkg_fs);

      // assuming field-rank 2
      atlas::Field stc_field = bkg_fs["stc"];
      int frank = stc_field.rank();
      if (frank > 2){
          oops::Log::error() << "Field stc rank "<< frank<<" but expected 2." << std::endl;
          throw eckit::BadValue("Erroneous field rank", Here());
      }

      // check required fields exist in the fieldset
      if (!bkg_fs.has("sheleg") ) {
          oops::Log::error() << "Missing required fields SWE in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_swe = atlas::array::make_view<double, 2>(bkg_fs["sheleg"]);
      if (!bkg_fs.has("vtype") ) {
          oops::Log::error() << "Missing required fields vtype in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_vtype = atlas::array::make_view<double, 2>(bkg_fs["vtype"]);
      if (!bkg_fs.has("stype") ) {
          oops::Log::error() << "Missing required fields stype in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_stype = atlas::array::make_view<double, 2>(bkg_fs["stype"]);

      if (!bkg_fs.has("stc") ) {
          oops::Log::error() << "Missing required fields stc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_stc = atlas::array::make_view<double, 2>(bkg_fs["stc"]);
      if (!bkg_fs.has("slc") ) {
          oops::Log::error() << "Missing required fields slc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_slc = atlas::array::make_view<double, 2>(bkg_fs["slc"]);
      if (!bkg_fs.has("soilMoistureVolumetric") ) {  // smc") ) {
          oops::Log::error() << "Missing required fields smc in state FieldSet. Aborting." << std::endl;
          throw eckit::BadValue("Missing required fields in state FieldSet", Here());
      }
      auto bkgv_smc = atlas::array::make_view<double, 2>(bkg_fs["soilMoistureVolumetric"]);
      // oops::Log::info() << "Finished reading background state" << std::endl;

      // vtype and stype to integers
      int len_land_vec = bkg_fs["sheleg"].shape(0);
      int lsoil = bkg_fs["stc"].shape(1);
      if (lsoil != lsoilc) {
	 oops::Log::error() << "lsoil " << lsoil << " must be equal to " << lsoilc << std::endl;
	 throw eckit::BadValue("The expected number of soil layers is 4 ", Here());
      }
      // state vectors. TODO: do this in the fort-cpp interface
      auto bkg_swe = viewToVector1D(bkgv_swe);
      auto bkg_vtype = viewToVector1D(bkgv_vtype);
      auto bkg_stype = viewToVector1D(bkgv_stype); 
      auto bkg_stc = viewToVector2D(bkgv_stc);
      auto bk_bkg_stc = viewToVector2D(bkgv_stc);
      auto bkg_slc = viewToVector2D(bkgv_slc);
      auto bkg_smc = viewToVector2D(bkgv_smc); 

      // Read increment
      // oops::Log::info() << "Reading increment" << std::endl;
      const eckit::LocalConfiguration incParams(fullConfig, "increment");
      oops::Variables incVars(incParams, "variables");
      // ijedi::Increment<ijedi::Traits> dx(geom_, incVars, xx.validTime());
      fv3jedi::Increment dx(geom_, incVars, xx.validTime());
      dx.read(incParams);
      oops::Log::test() << "Increment: " << dx << std::endl;

      bool upd_stc = false, upd_slc = false;
      atlas::FieldSet inc_fs;
      dx.toFieldSet(inc_fs);
      int lsoil_incr = incParams.getInt("lsoil_incr");
      std::vector<std::vector<double>> stc_inc(len_land_vec, std::vector<double>(lsoil_incr, 0.0f));
      if (inc_fs.has("stc")) {
        auto stcv_inc = atlas::array::make_view<double, 2>(inc_fs["stc"]);
        stc_inc = viewToVector2D(stcv_inc);
        upd_stc = true;
        oops::Log::info() << "Updating stc" << std::endl;
      }
      
      std::vector<std::vector<double>> slc_inc(len_land_vec, std::vector<double>(lsoil_incr, 0.0f));
      if (inc_fs.has("slc")) {
        auto slcv_inc = atlas::array::make_view<double, 2>(inc_fs["slc"]);
        slc_inc = viewToVector2D(slcv_inc);
        upd_slc = true;
        oops::Log::info() << "Updating slc" << std::endl;
      }
      // oops::Log::info() << "Done reading increment. len_land_vec: " << len_land_vec << std::endl;
 
      // read/construct mask for landice and snow tiles
      // std::vector<int> mask_landice(geom_.nlevsfc(), 0);
      std::vector<int> soil_mask(len_land_vec, 0);
      // TODO: check if landfrac and icefrac are relevant for mask
      SoilIncrements::calculateLandIncrementMask(bkg_swe, bkg_vtype, bkg_stype, 
                              len_land_vec, veg_type_landice, soil_mask);
      // zero out increments for mask not equal to 1
      for (int i = 0; i < len_land_vec; ++i) {
          if (soil_mask[i] != 1) {
              for (int j = 0; j < lsoil_incr; ++j) {
                  stc_inc[i][j] = 0.0;
                  slc_inc[i][j] = 0.0;
              }
          }
      }

      // Add increment to state
      // TODO: see if we can sync xx and dx configs, zero-out dx, then xx += dx;
      bool print_summary = false, print_debug = false;
      fullConfig.get("print_summary", print_summary);
      fullConfig.get("print_debug", print_debug);
      std::vector<int> stc_updated(len_land_vec, 0);
      std::vector<int> slc_updated(len_land_vec, 0);
      SoilIncrements::addIncrementSoil(
          myrank, lsoil, lsoil_incr, len_land_vec, soil_mask, 
          upd_stc, upd_slc, print_summary, print_debug,
          bkg_stc, bkg_slc, bkg_smc, stc_inc, slc_inc,       
          stc_updated, slc_updated  
      );
       
      // post-increment adjustments to ensure consistency b/n soil T and soil M
      SoilIncrements::applyLandDAadjustmentsSoil(
          lsoil_incr, isot, ivegsrc, len_land_vec, lsoil,
          istype, soil_mask, 
          bk_bkg_stc, bkg_stc, bkg_smc, bkg_slc,
          stc_updated, slc_updated, zsoil,
          upd_stc, upd_slc, myrank, print_summary, print_debug
      );

      // update state
      for (size_t i = 0; i < bkg_stc.size(); ++i) {
        for (size_t j = 0; j < bkg_stc[i].size(); ++j) {
          bkgv_stc(i, j) = static_cast<double>(bkg_stc[i][j]);
          bkgv_slc(i, j) = static_cast<double>(bkg_slc[i][j]);
	  bkgv_smc(i, j) = static_cast<double>(bkg_smc[i][j]);
        }
      }

      xx.fromFieldSet(bkg_fs);      
      oops::Log::test() << "Updated State: " << xx << std::endl;

      // Write updated state to file
      xx.write(eckit::LocalConfiguration(fullConfig, "output state"));

      return 0;
    }

   private:
      static constexpr std::array<float, 4> zsoil = { -0.1, -0.4, -1.0, -2.0 };
      static constexpr int veg_type_landice = 15;
      static constexpr int lsoilc = 4;     // zsoil is hard-coded for 4 layers
      static constexpr int ivegsrc = 1;   // The NOAHMP LSM expects that the ivegsrc physics parameter is 1
      static constexpr int isot = 1;      // Noahmp expects 1
      // hard coded defaults--unlikely to change
      bool frac_grid = true;
      float fice_threshold = 0.0;
      float lfrac_threshold = 0.0001;
      
      std::vector<double> viewToVector1D(const atlas::array::ArrayView<double, 2>& view) const {
        std::vector<double> vector1D(view.shape(0),0.0);
        for (size_t i = 0; i < view.shape(0); ++i) {
            vector1D[i] = view(i, 0); 
        }
        return vector1D;
      }  

       std::vector<std::vector<double>> viewToVector2D(const atlas::array::ArrayView<double, 2>& view) const {
	std::vector<std::vector<double>> vector2D(view.shape(0), std::vector<double>(view.shape(1),0.0));
        for (size_t i = 0; i < view.shape(0); ++i) {
          for (size_t j = 0; j < view.shape(1); ++j) {
            vector2D[i][j] = view(i, j);
          }
        }
        return vector2D;
      }

      // -----------------------------------------------------------------------------
      std::string appname() const override {
        return "landincr::AddLandIncrement";
      }
  };
}  // namespace land-apply_jedi_incr
