#ifndef __ACNetworkBlock
#define __ACNetworkBlock

#include "Block.h"
#include "LinearFunction.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "NetworkBlock.h"
#include "FRealObjective.h"

#include "wrapper_flow_cutter.h"

namespace SMSpp_di_unipi_it
{
class ACNetworkBlock : public NetworkBlock
{

 public:
 enum line_type
 {
  kNone = 0 ,  ///< no line
  kAC ,        ///< AC lines
  kHVDC ,      ///< HVDC lines
  kAC_HVDC     ///< AC and HVDC lines
 };

 class ACNetworkData : public NetworkBlock::NetworkData
 {
  public:
  ACNetworkData( void ) : NetworkBlock::NetworkData() {}
  ACNetworkData( NetworkData * ac_network_data ) {}
  virtual ~ACNetworkData() = default;
  virtual void deserialize( const netCDF::NcGroup & group ) override;
  Index get_number_lines( void ) const { return( f_number_lines ); }
  Index get_reference_node() const { return( f_reference_node ); }
  const std::vector< Index > & get_start_line( void ) const { return( v_start_line ); }
  const std::vector< Index > & get_end_line( void ) const { return( v_end_line ); }
  const std::vector< double > & get_min_power_flow( void ) const { return( v_min_power_flow ); }

  double get_min_power_flow( Index line ) const {
   assert( line < get_number_lines() );
   if( v_min_power_flow.empty() )
    return( 0 );
   return( v_min_power_flow[ line ] );
  }

  const std::vector< double > & get_max_power_flow( void ) const { return( v_max_power_flow ); }

  double get_max_power_flow( Index line ) const {
   assert( line < get_number_lines() );
   if( v_max_power_flow.empty() )
    return( 0 );
   return( v_max_power_flow[ line ] );
  }

  const std::vector< double > & get_susceptance( void ) const { return( v_susceptance ); }
  const std::vector< double > & get_network_cost( void ) const { return( v_network_cost ); }

  line_type get_lines_type( void ) const {
   if( get_number_lines() == 0 )
    return( kNone );
   if( std::all_of( v_susceptance.cbegin() , v_susceptance.cend() ,
                    []( double s ) { return( s == 0.0 ); } ) )
    return( kHVDC );
   if( std::all_of( v_susceptance.cbegin() , v_susceptance.cend() ,
                    []( double s ) { return( s != 0.0 ); } ) )
    return( kAC );
   return( kAC_HVDC );
  }

  const std::vector< std::string > & get_node_names( void ) const { return( v_node_names ); }
  const std::vector< std::string > & get_line_names( void ) const { return( v_line_names ); }

  void run_tree_decomposition(){
    /* in prevision */
    std::vector<int> start_line(std::begin(v_start_line), std::end(v_start_line));
    std::vector<int> end_line(std::begin(v_end_line), std::end(v_end_line));
    run_flow_cutter(
            get_number_nodes(), 
            start_line, 
            end_line, 
            v_tree_bags,
            v_start_tree_edges,
            v_end_tree_edges,
            "min_shortcut");
  }

  virtual void serialize( netCDF::NcGroup & group ) const override;

  protected:

  Index f_number_lines{};
  Index f_reference_node;    ///< reference node (used in the PTDF matrix)
  std::vector< Index > v_start_line;
  std::vector< Index > v_end_line;
  std::vector< double > v_susceptance;
  std::vector< double > v_min_power_flow;
  std::vector< double > v_max_power_flow;
  std::vector< double > v_network_cost;
  std::vector< std::string > v_node_names;  ///< Node names
  std::vector< std::string > v_line_names;  ///< Line names

  std::vector< std::vector<int> > v_tree_bags;
  std::vector< int > v_start_tree_edges;
  std::vector< int > v_end_tree_edges;

  private:

  SMSpp_insert_in_factory_h;

 };  // end( class( ACNetworkData ) )

 explicit ACNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_NetworkData( nullptr ) {}

 virtual ~ACNetworkBlock() override;

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

 void generate_objective( Configuration * objc = nullptr ) override;
 
 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_number_nodes() );
 }

 Index get_number_lines( void ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_number_lines() );
 }

 double get_kappa( Index line ) const {
  if( v_kappa.empty() )
   return( 1 );
  assert( line < v_kappa.size() );
  return( v_kappa[ line ] );
 }

 double get_min_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_min_power_flow( line ) );
 }

 double get_max_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_max_power_flow( line ) );
 }

 NetworkData * get_NetworkData( void ) const override {
  return( f_NetworkData );
 }

 const double * get_active_demand( Index interval = 0 ) const override {
  if( v_ActiveDemand.empty() )
   return( nullptr );
  return( &( v_ActiveDemand.front() ) );
 }

 const std::vector< ColVariable > & get_power_flow( void ) const {
  return( v_power_flow );
 }

 const std::vector< ColVariable > & get_auxiliary_variable( void ) const {
  return( v_auxiliary_variable );
 }

 const std::vector< FRowConstraint > &
 get_power_flow_limit_constraints( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error( "ACNetworkBlock::get_power_flow_limit_constraints:"
                            " ACNetworkData has not been set." ) );

  switch( f_NetworkData->get_lines_type() ) {
   case( kAC ):
    return( v_AC_power_flow_limit_const );
   case( kAC_HVDC ):
   default:
    return( v_AC_HVDC_power_flow_limit_const );
  }
 }

 const std::vector< BoxConstraint > &
 get_power_flow_limit_HVDC_bounds( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error( "ACNetworkBlock::get_power_flow_limit_HVDC_bounds:"
                            " ACNetworkData has not been set." ) );
  return( v_HVDC_power_flow_limit_const );
 }

 void set_NetworkData( NetworkBlock::NetworkData * nd = nullptr ) override {
  // if there was a previous ACNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = static_cast< ACNetworkData * >( nd );
  f_local_NetworkData = false;
 }

 void set_ActiveDemand(
  const std::vector< std::vector< double > > & v ) override {
  if( v_ActiveDemand.empty() )
   v_ActiveDemand = v[ 0 ];
 }

 void deserialize( const netCDF::NcGroup & group ) override;

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "ACNetworkBlock::load() not implemented yet" ) );
 }

 void serialize( netCDF::NcGroup & group ) const override;

 void run_tree_decomposition();

 void set_kappa( MF_dbl_it values ,
                 Subset && subset ,
                 const bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

 void change_power_flow_limit_constraints( 
                const std::vector<Index>& modified_lines, 
                c_ModParam issueAMod);

 void change_relax_abs_constraints(
                const std::vector<Index>& modified_lines, 
                c_ModParam issueAMod);

 void change_DC_power_flow_injection_constraints(
                const std::vector<Index>& modified_nodes, 
                c_ModParam issueAMod);

 void set_active_demand( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 void set_active_demand( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 protected:

 ACNetworkData * f_NetworkData;
 std::vector< double > v_ActiveDemand;
 std::vector< double > v_kappa;
 std::vector< ColVariable > v_power_flow;
 std::vector< ColVariable > v_auxiliary_variable;
 std::vector< FRowConstraint > v_AC_power_flow_limit_const;
 std::vector< FRowConstraint > v_AC_HVDC_power_flow_limit_const;
 std::vector< FRowConstraint > v_power_flow_injection_const;
 boost::multi_array< FRowConstraint , 2 > v_power_flow_relax_abs;
 std::vector< BoxConstraint > v_HVDC_power_flow_limit_const;
 std::vector< BoxConstraint > node_injection_bounds_const;
 FRealObjective objective;

 private:

 SMSpp_insert_in_factory_h;

 static void static_initialization( void ) {

  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "ACNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "ACNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );
 }

};  // end( class( ACNetworkBlock ) )


class ACNetworkBlockMod : public NetworkBlockMod
{

 public:

 /// public enum for the types of ACNetworkBlockMod
 enum ACNetB_mod_type
 {
  eSetKappa = eNetBModLastParam ,  ///< set the kappa constants
  eACNetBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of ACNetworkBlockMod. */
 };

 /// constructor, takes the ACNetworkBlock and the type
 ACNetworkBlockMod( ACNetworkBlock * const fblock , const int type )
  : NetworkBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~ACNetworkBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the ACNetworkBlockMod
 void print( std::ostream & output ) const override {
  output << "ACNetworkBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set active demand values ";
  }
 }

 ACNetworkBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( ACNetworkBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ACNetworkBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ACNetworkBlockMod for "ranged" modifications
class ACNetworkBlockRngdMod : public ACNetworkBlockMod
{

 public:

 /// constructor: takes the ACNetworkBlock, the type, and the range
 ACNetworkBlockRngdMod( ACNetworkBlock * const fblock , const int type ,
                        Block::Range rng )
  : ACNetworkBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~ACNetworkBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the ACNetworkBlockRngdMod
 void print( std::ostream & output ) const override {
  ACNetworkBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( ACNetworkBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ACNetworkBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ACNetworkBlockMod for "subset" modifications
class ACNetworkBlockSbstMod : public ACNetworkBlockMod
{

 public:

 /// constructor: takes the ACNetworkBlock, the type, and the subset
 ACNetworkBlockSbstMod( ACNetworkBlock * const fblock , const int type ,
                        Block::Subset && nms )
  : ACNetworkBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ACNetworkBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the ACNetworkBlockSbstMod
 void print( std::ostream & output ) const override {
  ACNetworkBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( ACNetworkBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ACNetworkBlock.h included */