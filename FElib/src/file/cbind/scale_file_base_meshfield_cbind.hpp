#pragma once

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "../../mesh/cbind/scale_mesh_linedom1d_cbind.hpp"
#include "../../mesh/cbind/scale_mesh_base1d_cbind.hpp"

#include "../../common/cbind/scale_common_cbind_util.hpp"
using namespace scale::FElib;

extern "C" {
    void* CFILE_base_meshfield_Init(int var_num, 
        void* mesh1D_ptr, void* mesh2D_ptr, void* meshCubedSphere2D_ptr, void* mesh3D_ptr, void* meshCubedSphere3D_ptr,
        bool force_uniform_grid);
    void CFILE_base_meshfield_Final(void* ptr);

    void CFILE_base_meshfield_Open(void* ptr, const char* basename, int myrank);
    void CFILE_base_meshfield_Create(void* ptr, const char* basename, const char* title, const char* dtype, 
        bool* fileexisited, int myrank, const char* tunits, const char* calendar );


    void CFILE_base_meshfield_Close(void* ptr);

    void CFILE_base_meshfield_Def_var1(void* ptr, void* field, const char* desc, int vid, int dim_type_id, const char* datatype, const char* standard_name, double timeinv, int nsteps);
    void CFILE_base_meshfield_End_def(void* ptr);

    void CFILE_base_meshfield_Write_var1d(void* ptr, int vid, void* field1d, double sec_str, double sec_end);
}
class File_base_meshfield{
public:
    File_base_meshfield(int varnum, const MeshBase1D& mesh1D, bool force_uniform_grid=false)
        : handle_(cbind::init_handle(
              &CFILE_base_meshfield_Init, &CFILE_base_meshfield_Final, "Failed to initialize CFile_base_meshfield",
              varnum, mesh1D.get_Handle().get(), nullptr, nullptr, nullptr, nullptr, force_uniform_grid) ){}

    ~File_base_meshfield() = default;

    File_base_meshfield(const File_base_meshfield&) = delete;
    File_base_meshfield& operator=(const File_base_meshfield&) = delete;
    File_base_meshfield(File_base_meshfield&& other) noexcept = default;
    File_base_meshfield& operator=(File_base_meshfield&& other) = default;

    void open(const std::string& basename, int myrank=-1){ CFILE_base_meshfield_Open(handle_.get(), basename.c_str(), myrank); }

    void create(const std::string& basename, const std::string& title, const std::string& dtype, bool& fileexisted, int myrank, const std::string& tunits, const std::string& calendar){
        CFILE_base_meshfield_Create( handle_.get(), basename.c_str(), title.c_str(), dtype.c_str(), &fileexisted, myrank, tunits.c_str(), calendar.c_str());
    }
    void create(const std::string& basename, const std::string& title, const std::string& dtype, bool& fileexisted, int myrank=-1){
        CFILE_base_meshfield_Create( handle_.get(), basename.c_str(), title.c_str(), dtype.c_str(), &fileexisted, myrank, "seconds", "");
    }
    
    void close(){ CFILE_base_meshfield_Close(handle_.get()); }

    void def_var(const MeshFieldBase& field, const std::string& desc, int vid, int dim_type_id, const std::string& datatype, const std::string& standard_name, double timeinv, int nsteps){ 
        CFILE_base_meshfield_Def_var1(handle_.get(), field.get_Handle().get(), desc.c_str(), vid, dim_type_id, datatype.c_str(), standard_name.c_str(), timeinv, nsteps); }

    void end_def(){ CFILE_base_meshfield_End_def(handle_.get()); }

    void write_var(int vid, const MeshField1D& field, double sec_str, double sec_end){ CFILE_base_meshfield_Write_var1d(handle_.get(), vid, field.get_Handle().get(), sec_str, sec_end); }

private:
    cbind::Handle handle_;
};
