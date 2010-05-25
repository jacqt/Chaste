"""Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
"""

"""
Defines the Protocol class, which encapsulates the input & output of a
simulation protocol.
"""

import pycml
from pycml import *
import validator

class ProtocolError(ValueError):
    """Error thrown if a Protocol instance is invalid."""
    pass

class Protocol(object):
    """A class representing a simulation protocol.

     * When a protocol is initialised, it should be passed a cellml_model, and
       alter the model equations to reflect the protocol: remove ODEs for
       state variables that are set by the protocol, replace equations for
       protocol inputs with forms specified by the protocol.  Then redo the
       topological sort to ensure the equations are ordered correctly?
     * The main translator classes should then not need to do anything special
       when using a protocol, except for considering the case where there are
       no ODEs.
     * Classes need a ComputeDerivedQuantities method, which takes an (optional)
       state variable vector, and time, and computes all quantities in the model
       which are non-state-var outputs.  The protocol system, possibly via
       model annotations, will need to be able to specify what goes in here.  We
       may need the ComputedVariable functionality in OdeSystemInformation to
       name entries in the result vector.  The pe:keep behaviour for computed
       vars could change to store variables in this vector?  Although we don't
       necessarily want all ionic currents in here - the main reason they get
       annotated is so that we can generate GetIIonic.  So I think we need
       separate annotations for specifying parameters and computed vars (which
       both imply pe:keep).
     * GetIIonic may need a time parameter, since protocol inputs may depend on
       time, and the generated code wouldn't compile otherwise.
    
    A fully initialised protocol contains the following attributes:
     * inputs - a list of protocol inputs.  These may be cellml_variable instances,
       to (re)define a variable in the model, or mathml_apply instances, to add or
       modify an equation.  Once modify_model has been called, these objects will
       also exist in the model.
    """
    def __init__(self, model, multi_stage=False):
        """Create a new protocol.
        
        Eventually this will have arguments to parse a protocol definition file
        and create the protocol.  For now, however, the only option is to specify
        multi_stage as True and set up the internal data structures yourself.
        Then call self.modify_model.
        """
        self._protocol_component = None
        self.model = model
        self.inputs = []
        if not multi_stage:
            raise NotImplemented

    def modify_model(self):
        """Actually apply protocol modifications to the model.
        
        Prior to this being called, all variable references within self.inputs must
        use full names, i.e. 'component_name,variable_name'.  Variables without a
        component part will be placed into a new 'protocol' component.
        
        This method will add the items from self.inputs into the model, replacing
        variables with the same name, and equations that assign to the same variable.
        This may involve changing a variable's type from State to Computed, or vice
        versa.
        
        After the call, all names and name references will be 'local'.  Connections
        will be created between components as needed, and units definitions added, to
        ensure a valid model.  In order for this to work, if a variable has units that
        do not already exist in the model, the object *must* have an attribute
        _cml_units referring to a suitable cellml_units instance.
        """
        for input in self.inputs:
            if isinstance(input, cellml_variable):
                self._add_variable_to_model(input)
            elif isinstance(input, mathml_apply):
                self._add_maths_to_model(input)
            else:
                raise ProtocolError("Unexpected input object type '" + str(type(input))
                                    + "': " + str(input))
        self._clear_model_caches()
        self._fix_model_connections()
        self._reanalyse_model()
        
    def _reanalyse_model(self):
        """Re-do the model validation steps needed for further processing of the model.
        
        Checks connections, etc. and builds up the dependency graph again, then performs
        a topological sort.
        """
        # We want to see any errors
        logging_info = validator.CellMLValidator.setup_logging(show_errors=True, show_warnings=False)
        # Re-run validation & analysis
        self.model._check_variable_mappings()
        if not self.model._cml_validation_errors:
            assignment_exprs = self.model.search_for_assignments()
            self.model._check_assigned_vars(assignment_exprs)
        if not self.model._cml_validation_errors:
            self.model._classify_variables(assignment_exprs)
            self.model._order_variables(assignment_exprs)
        if self.model._cml_validation_errors:
            raise ProtocolError("Applying protocol created an invalid model.")
        # Clear up logging
        validator.CellMLValidator.cleanup_logging(logging_info)
    
    def _fix_model_connections(self):
        """Ensure the modified model has all the necessary connections between variables.
        
        Check mathematics for ci elements that refer to variables not defined in that
        component.  These must refer to the variable by its full name (i.e. 'cname,vname').
        These variables will be renamed to use local names by this method, which will
        also create local variables mapped to the relevant source variable if needed.
        
        This needs to take account of the fact that a variable in one nested component
        may need to be connected to a variable in another nested component, and so create
        variables in the parent components to connect the whole thing up.
        """
        for expr in self.model.search_for_assignments():
            for ci_elt in self._find_ci_elts(expr):
                vname = unicode(ci_elt)
                if u',' in vname:
                    cname, vname = self._split_name(vname)
                    comp = expr.component
                    if comp.name != cname:
                        raise NotImplemented
                    else:
                        # Just rename to be local
                        ci_elt._rename(vname)
    
    def _clear_model_caches(self):
        """
        Clear cached links in the model, since we'll need to recompute many of them
        once we've finished modifying it.  Also clears dependency information.
        """
        for comp in getattr(self.model, u'component', []):
            for math in getattr(comp, u'math', []):
                math._unset_cached_links()
        for var in self.model.get_all_variables():
            var.clear_dependency_info()
        assignment_exprs = self.model.search_for_assignments()
        for expr in assignment_exprs:
            expr.clear_dependency_info()
    
    def _get_protocol_component(self):
        """Get the protocol component in the model, creating it if necessary.
        
        New variables created just for use by the simulation protocol get put into a
        new 'protocol' component.  If a component with that name already exists,
        underscores will be added to the component name to make it unique.
        """
        if self._protocol_component is None:
            cname = u'protocol'
            while True:
                try:
                    comp = self.model.get_component_by_name(cname)
                    cname += u'_'
                except KeyError:
                    # Component with this name doesn't exist
                    break
            # Create the component
            comp = cellml_component.create_new(self.model, cname)
            self.model._add_component(comp)
            self._protocol_component = comp
        return self._protocol_component
    
    def _split_name(self, full_name):
        """Split a full name into cname,vname, creating the component if needed.
        
        If the full_name doesn't contain a component part, the 'protocol' component
        will be used.
        """
        parts = full_name.split(',')
        if len(parts) == 2:
            cname, vname = parts
        elif len(parts) == 1:
            cname = self._get_protocol_component().name
            vname = full_name
        else:
            raise ValueError("Invalid variable name: " + full_name)
        return cname, vname
        
    def _rename_local_variables(self, expr):
        """
        Change local variable references in the given expression to refer
        explicitly to the protocol component.
        """
        for ci_elt in self._find_ci_elts(expr):
            vname = unicode(ci_elt)
            if u',' not in vname:
                # Do the rename
                cname = self._get_protocol_component().name
                full_name = cname + u',' + vname
                expr._rename(full_name)
    
    def _find_ci_elts(self, expr):
        """Get an iterator over all ci elements on the descendent-or-self axis of the given element."""
        if isinstance(expr, mathml_ci):
            yield expr
        elif hasattr(expr, 'xml_children'):
            # Recurse
            for child in expr.xml_children:
                for ci_elt in self._find_ci_elts(child):
                    yield ci_elt
        
    def _add_variable_to_model(self, var):
        """Add or replace a variable in our model.
        
        We don't really do any checking for this case - just add the variable.
        This means that some 'possible' changes don't actually make sense, for
        instance giving a computed variable an initial value will trigger a later
        validation error, unless its definition is also changed to an ODE.
        (To change a variable to a constant, you need to replace its definition
        with a constant expression.)
        """
        cname, vname = self._split_name(var.name)
        comp = self.model.get_component_by_name(cname)
        orig_var = comp.get_variable_by_name(vname)
        if orig_var:
            # We're replacing a variable
            comp._del_variable(orig_var)
        var.name = vname
        comp._add_variable(var)

    def _add_maths_to_model(self, expr):
        """Add or replace an equation in the model.
        
        This case is more complex than variables, since we may need to change
        the type of variable assigned to, depending on the expression.
        
        Note: variable references within ci elements in the given expression
        should use full names (i.e. cname,vname).  Any local names will be
        assumed to refer to variables in the protocol component, and modified
        by self._rename_local_variables.  Later, self._fix_model_connections
        will change all references to use local names.
        """
        assert isinstance(expr, mathml_apply)
        assert expr.operator().localName == u'eq', 'Expression is not an assignment'
        # Figure out what's on the LHS of the assignment
        lhs = expr.operands().next()
        if lhs.localName == u'ci':
            # Straight assignment to variable
            cname, vname = self._split_name(unicode(lhs))
            assigned_var = self.model.get_variable_by_name(cname, vname)
            self._remove_existing_definition(assigned_var, False)
            self._add_expr_to_comp(cname, expr)
        else:
            # This had better be an ODE
            assert lhs.localName == u'apply', 'Expression is not a straight assignment or ODE'
            assert lhs.operator().localName == u'diff', 'Expression is not a straight assignment or ODE'
            dep_var = lhs.operands().next()
            assert dep_var.localName == u'ci', 'ODE is malformed'
            cname, dep_var_name = self._split_name(unicode(dep_var))
            dep_var = self.model.get_variable_by_name(cname, dep_var_name)
            self._remove_existing_definition(dep_var, True)
            self._add_expr_to_comp(cname, expr)

    def _remove_existing_definition(self, var, keep_initial_value):
        """Remove any existing definition (as an equation) of the given variable.
        
        If keep_initial_value is False, then also remove any initial_value attribute.
        
        If the variable is Mapped, throw a ProtocolError.
        """
        if var.get_type() == VarTypes.Mapped:
            raise ProtocolError("Cannot add new mathematics defining a mapped variable - change the definition of its source instead")
        if not keep_initial_value and hasattr(var, u'initial_value'):
            del var.initial_value
        # Note: if this is a variable added by the protocol, then it shouldn't have
        # any dependencies set up yet, so this is a no-op.
        for dep in var._get_all_expr_dependencies():
            assert isinstance(dep, mathml_apply)
            dep.xml_parent.xml_remove_child(dep)
            dep.xml_parent = None # Not done by Amara...
    
    def _add_expr_to_comp(self, cname, expr):
        """Add an expression to the mathematics in the given component."""
        comp = self.model.get_component_by_name(cname)
        if not hasattr(comp, u'math'):
            # Create the math element
            math = comp.xml_create_element(u'math', NSS[u'm'])
            comp.xml_append(math)
        # Append this expression
        comp.math.xml_append(expr)
