#ifndef CARDIACELECTROMECHANICSPROBLEM_HPP_
#define CARDIACELECTROMECHANICSPROBLEM_HPP_

#include "MooneyRivlinMaterialLaw.hpp"
#include "CardiacMechanicsAssembler.cpp"
#include "FiniteElasticityTools.hpp"


/**
 *  Solve a full cardiac electro-mechanics problem in 2d or 3d.
 * 
 *  See documentation for AbstractCardiacElectroMechanicsProblem
 */
template<unsigned DIM>
class CardiacElectroMechanicsProblem : public AbstractCardiacElectroMechanicsProblem<DIM>
{
private:
    unsigned mNumElementsPerDimInElectricsMesh;
    unsigned mNumElementsPerDimInElectricsMesh;

public:
    /** 
     *  Constructor
     *  @param pCellFactory cell factory for creating cells (see Monodomain tests)
     *  @endTime end time of the simulation. Start time is assumed to be 0.0
     *  @timeStep time step for the electrics (and currently the mechanics too)
     *  @useExplicit Whether to use an explicit or implicit mechanics solver
     *  @outputDirectory. Output directory. Omit if no output is required.
     * 
     *  See documentation for AbstractCardiacElectroMechanicsProblem
     */
    CardiacElectroMechanicsProblem(AbstractCardiacCellFactory<DIM>* pCellFactory,
                                   double endTime,
                                   double timeStep,
                                   bool useExplicitMethod,
                                   unsigned numElementsPerDimInElectricsMesh,
                                   unsigned numElementsPerDimInMechanicsMesh,
                                   std::string outputDirectory = "")
        :  AbstractCardiacElectroMechanicsProblem<DIM>(pCellFactory,
                                                       endTime,
                                                       timeStep,
                                                       useExplicitMethod,
                                                       outputDirectory)
    {
        assert(numElementsPerDimInElectricsMesh > 8);
        numElementsPerDimInElectricsMesh
        numElementsPerDimInElectricsMesh = numElementsPerDimInElectricsMesh;
        numElementsPerDimInElectricsMesh = numElementsPerDimInElectricsMesh
        
    }
    

    void ConstructMeshes()
    {        
        // create electrics mesh
        this->mpElectricsMesh = new ConformingTetrahedralMesh<DIM,DIM>();

        unsigned num_elem = 40;
        this->mpElectricsMesh->ConstructRectangularMesh(num_elem,num_elem);
        this->mpElectricsMesh->Scale(1.0/num_elem,1.0/num_elem);

        // create mechanics mesh
        this->mpMechanicsMesh = new Triangulation<DIM>();
        GridGenerator::hyper_cube(*(this->mpMechanicsMesh), 0.0, 1.0);
        this->mpMechanicsMesh->refine_global(3);
        
        std::cout << "numnodes = " << this->mpElectricsMesh->GetNumNodes() << ", " << this->mpMechanicsMesh->n_vertices() << "\n";
        
        //assert(this->mpMechanicsMesh->n_vertices()<=this->mpElectricsMesh->GetNumNodes());
    }

    
    void ConstructMechanicsAssembler(std::string mechanicsOutputDir)
    {
        Point<DIM> zero;
        FiniteElasticityTools<DIM>::FixFacesContainingPoint(*(this->mpMechanicsMesh), zero);
        
        MooneyRivlinMaterialLaw<DIM>* p_law = new MooneyRivlinMaterialLaw<DIM>(2.0);
        if(this->mUseExplicitMethod)
        {
            this->mpCardiacMechAssembler = new CardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh,mechanicsOutputDir,p_law);
        }
        else
        {
            assert(0); // not done yet..
            //this->mpCardiacMechAssembler = new ImplicitCardiacMechanicsAssembler<DIM>(this->mpMechanicsMesh);
        }
    }
};

#endif /*CARDIACELECTROMECHANICSPROBLEM_HPP_*/
