#include "solver.hpp"
#include "init.h"

Solver::Solver(Mesh& mesh, Parameter& param)
    : mesh_(mesh), param_(param)
{
    fieldCreator_ = new FieldCreator(mesh_, param_);
    gl_ = new Gauss_Lobatto(param_.numGLP);
    baseFunction_ = new BaseFunction(param_.dimPk, param_.numGLP, *gl_, mesh_.dx(), mesh_.dy());

    output_ = new Output(mesh_);
}

Solver::~Solver()
{
    if (fieldCreator_) delete fieldCreator_;
    if (gl_) delete gl_;
    if (baseFunction_) delete baseFunction_;
    if (output_) delete output_;
}

void Solver::init()
{
    // create fields
    createFields();

    // set initial condition
    setInitialCondition();
}

void Solver::createFields()
{
    fields_.clear();

    // 守恒变量
    fields_.emplace("rho", fieldCreator_->createField());
    fields_.emplace("rhou", fieldCreator_->createField());
    fields_.emplace("rhov", fieldCreator_->createField());
    fields_.emplace("E", fieldCreator_->createField());

    // 原始变量
    fields_.emplace("u", fieldCreator_->createField());
    fields_.emplace("v", fieldCreator_->createField());
    fields_.emplace("p", fieldCreator_->createField());

    // Field& u = fields.at("u");
    // auto [nx, ny, dimPk] = u.dataShape();
    // std::cout << "nx = " << nx << ", ny = " << ny << ", dimPk = " << dimPk << std::endl;
}

void Solver::setInitialCondition()
{

    size_t nx = mesh_.nx();  size_t ny = mesh_.ny();
    scalar dx = mesh_.dx();  scalar dy = mesh_.dy();
    scalar dx1 = dx * 0.5;   scalar dy1 = dy * 0.5;

    const std::vector<scalar>& xc = mesh_.xc();
    const std::vector<scalar>& yc = mesh_.yc();

    const std::vector<scalar>& lambda = gl_->lambda();

    array_3& rho = fields_.at("rho").phi();
    array_3& u = fields_.at("u").phi();
    array_3& v = fields_.at("v").phi();
    array_3& p = fields_.at("p").phi();

    size_t numGLP = param_.numGLP;

    scalar x, y; // coord

    for(size_t i = 0; i < nx; i++)
    for(size_t j = 0; j < ny; j++)
    for(size_t k = 0; k < numGLP; k++)
    {
        x = xc[i] + lambda[k] * dx1;
        y = yc[j] + lambda[k] * dy1;

        rho[i][j][k] = rho0(x, y);
        u[i][j][k] = u0(x, y);
        v[i][j][k] = v0(x, y);
        p[i][j][k] = p0(x, y);
    } 
 
    output_->writeVtk("p", p);


}