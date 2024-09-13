#include "solver.hpp"
#include "init.h"

Solver::Solver(Mesh& mesh, Parameter& param)
    : mesh_(mesh), param_(param)
{
    fieldCreator_ = new FieldCreator(mesh_, param_);
    gl_ = new Gauss_Lobatto(param_.numGLP);
    baseFunction_ = new BaseFunction(param_.dimPk, param_.numGLP, *gl_, mesh_.dx(), mesh_.dy());
    phi2coff = new Phi2Coff(gl_->weight(), baseFunction_->phiGauss(), baseFunction_->mm());

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
    createFields();        // 创建场
    setInitialCondition(); // 设置初场
    computeL2forInitial(); // 计算初场L2投影，得到基函数系数初值
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

    label nx = mesh_.nx();  label ny = mesh_.ny();
    scalar dx = mesh_.dx();  scalar dy = mesh_.dy();
    scalar dx1 = dx * 0.5;   scalar dy1 = dy * 0.5;

    const std::vector<scalar>& xc = mesh_.xc();
    const std::vector<scalar>& yc = mesh_.yc();

    const std::vector<scalar>& lambda = gl_->lambda();

    array_4& rho = fields_.at("rho").phi();
    array_4& u = fields_.at("u").phi();
    array_4& v = fields_.at("v").phi();
    array_4& p = fields_.at("p").phi();

    label numGLP = param_.numGLP;

    scalar x, y; // coord

    for(label i = 0; i < nx; i++)
    for(label j = 0; j < ny; j++)
    for(label i1 = 0; i1 < numGLP; i1++)
    for(label j1 = 0; j1 < numGLP; j1++)
    {
        x = xc[i] + lambda[i1] * dx1;
        y = yc[j] + lambda[j1] * dy1;

        rho(i,j,i1,j1) = rho0(x, y);
        u(i,j,i1,j1) = u0(x, y);
        v(i,j,i1,j1) = v0(x, y);
        p(i,j,i1,j1) = p0(x, y);
    } 
 
    array_4& rhou = fields_.at("rhou").phi();
    array_4& rhov = fields_.at("rhov").phi();
    array_4& E = fields_.at("E").phi();

    rhou = rho * u;
    rhov = rho * v;
    E = p / (param_.gamma - 1) + 0.5 * rho * (u * u + v * v);


}

void Solver::computeL2forInitial()
{
    for(auto& field : fields_)
    {
        phi2coff->computeCoff(field.second);
    }

    array_3& p = fields_.at("E").data();
    output_->writeVtk("E", p);
}
