
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
This module abstracts the interface to RDF metadata about CellML models.
"""

import types

import RDF

import pycml

# Map from cellml_model instances to RDF.Model instances
_models = {}

# Base URI to use for models.  Unfortunately the RDF library won't let
# us use an empty URI, so we use a dummy URI then strip it out when
# serializing the RDF.
_base_uri = 'urn:chaste-pycml:dummy-rdf-base-uri'


# Allowed metadata names, more to come
# TODO #1209: Use a proper ontology!
METADATA_NAMES = frozenset(
    ['membrane_voltage', 'membrane_capacitance', 'membrane_E_R', 'membrane_stimulus_current', 'membrane_stimulus_current_duration',
     'membrane_stimulus_current_amplitude','membrane_stimulus_current_period','membrane_stimulus_current_offset', 'sodium_channel_current',
     'sodium_channel_conductance', 'sodium_channel_m_gate', 'sodium_channel_h_gate', 
     'potassium_channel_current', 'potassium_channel_conductance', 'potassium_channel_n_gate', 
     'leakage_current', 'sodium_channel_current_conductance','sodium_channel_current_h_gate',
     'sodium_channel_current_j_gate','sodium_channel_current_m_gate', 'temperature',
     'potassium_reversal_potential_sodium_permeability', 'inward_rectifier_potassium_current_conductance',
     'rapid_time_dependent_potassium_current_conductance', 'rapid_time_dependent_potassium_current_Xr1_gate',
     'rapid_time_dependent_potassium_current_Xr2_gate', 'slow_time_dependent_potassium_current_conductance',
     'slow_time_dependent_potassium_current_Xs_gate', 'fast_sodium_current_conductance',
     'fast_sodium_current_m_gate', 'fast_sodium_current_h_gate', 'fast_sodium_current_j_gate', 'sodium_background_current_conductance',
     'L_type_Ca_current_conductance', 'L_type_Ca_current_d_gate', 'L_type_Ca_current_f_gate', 'L_type_Ca_current_f2_gate',
     'L_type_Ca_current_fCass_gate', 'calcium_background_current_conductance', 'transient_outward_current_conductance',
     'transient_outward_current_s_gate', 'transient_outward_current_r_gate', 'sodium_potassium_pump_current_permeability',
     'sodium_calcium_exchanger_current_maximum', 'calcium_pump_current_conductance', 'potassium_pump_current_conductance',
     'calcium_dynamics_release_current_maximum', 'calcium_dynamics_leak_current_maximum', 'calcium_leak_current_conductance',
     'calcium_dynamics_uptake_current_maximum',
     'cytosolic_calcium_concentration'])



def _debug(*args):
    pycml.DEBUG('cellml-metadata', *args)


def _get_rdf_from_model(cellml_model):
    """Get the RDF model of the given CellML model.
    
    If this model is already in our map, return the existing RDF store.
    Otherwise, extract metadata from all RDF elements in the cellml_model,
    create a new RDF model from these, and delete the original elements.
    """
    if not cellml_model in _models:
        rdf_blocks = cellml_model.xml_xpath(u'//rdf:RDF')
        m = RDF.Model()
        p = RDF.Parser()
        for rdf_block in rdf_blocks:
            rdf_text = rdf_block.xml()
            p.parse_string_into_model(m, rdf_text, _base_uri)
            rdf_block.xml_parent.xml_remove_child(rdf_block)
        _models[cellml_model] = m
    return _models[cellml_model]

def remove_model(cellml_model):
    """The given model is being deleted / no longer needed."""
    if cellml_model in _models:
        del _models[cellml_model]
        _debug('Clearing RDF state for model', cellml_model.name)

def update_serialized_rdf(cellml_model):
    """Ensure the RDF serialized into the given CellML model is up-to-date.
    
    If we have done any metadata processing on the given model, will serialize
    our RDF store into the rdf:RDF element child of the model.
    """
    if cellml_model in _models:
        # Paranoia: ensure it doesn't already contain serialized RDF
        if hasattr(cellml_model, u'RDF'):
            pycml.LOG('cellml-metadata', logging.WARNING, 'Removing existing RDF in model.')
            old_rdf_block = cellml_model.RDF
            old_rdf_block.xml_parent.xml_remove_child(old_rdf_block)
        # Serialize the rdf model into cellml_model.RDF
        m = _models[cellml_model]
        rdf_text = m.to_string().replace(_base_uri, '')
        rdf_doc = pycml.amara.parse(rdf_text)
        cellml_model.xml_append(rdf_doc.RDF)
        # Remove the RDF model
        remove_model(cellml_model)

def create_rdf_node(node_content=None, fragment_id=None):
    """Create an RDF node.
    
    node_content, if given, must either be a tuple (qname, namespace_uri),
    or a string, in which case it is interpreted as a literal RDF node.
    
    Alternatively, fragment_id may be given to refer to a cmeta:id within the
    current model.
    
    If neither are given, a blank node is created.
    """
    if fragment_id:
        node = RDF.Node(uri_string=str(_base_uri+'#'+fragment_id))
    elif node_content:
        if type(node_content) == types.TupleType:
            qname, nsuri = node_content
            if nsuri[-1] not in ['#', '/']:
                nsuri = nsuri + '#'
            prefix, local_name = pycml.SplitQName(qname)
            node = RDF.Node(uri_string=str(nsuri+local_name))
        elif type(node_content) in types.StringTypes:
            node = RDF.Node(str(node_content))
        else:
            raise ValueError("Don't know how to make a node from " + str(node_content)
                             + " of type " + type(node_content))
    else:
        node = RDF.Node()
    return node

def replace_statement(cellml_model, source, property, target):
    """Add a statement to the model, avoiding duplicates.
    
    Any existing statements with the same source and property will first be
    removed.
    """
    _debug("replace_statement(", source, ",", property, ",", target, ")")
    rdf_model = _get_rdf_from_model(cellml_model)
    # Check for existing statements
    query = RDF.Statement(subject=source, predicate=property, object=None)
    for statement in rdf_model.find_statements(query):
        del rdf_model[statement]
    # Add the new statement
    statement = RDF.Statement(subject=source, predicate=property, object=target)
    rdf_model.append(statement)

def get_target(cellml_model, source, property):
    """Get the target of property from source.
    
    Returns None if no such target exists.  Throws if there is more than one
    match.
    
    If the target is a literal node, returns its string value.  Otherwise returns
    an RDF node.
    """
    rdf_model = _get_rdf_from_model(cellml_model)
    targets = list(rdf_model.targets(source, property))
    if len(targets) > 1:
        raise ValueError("Too many targets for source " + str(source) + " and property " + str(property))
    elif len(targets) == 1:
        target = targets[0]
    else:
        target = None
    if target and target.is_literal():
        target = str(target)
    _debug("get_target(", source, ",", property, ") -> ", "'" + str(target) + "'")
    return target

def find_variables(cellml_model, property, value=None):
    """Find variables in the cellml_model with the given property, and optionally value.
    
    property (and value if given) should be a suitable input for create_rdf_node.
    
    Will return a list of cellml_variable instances.
    """
    _debug("find_variables(", property, ",", value, ")")
    rdf_model = _get_rdf_from_model(cellml_model)
    property = create_rdf_node(property)
    if value:
        value = create_rdf_node(value)
    query = RDF.Statement(None, property, value)
    results = rdf_model.find_statements(query)
    vars = []
    for result in results:
        assert result.subject.is_resource(), "Non-resource annotated."
        uri = str(result.subject.uri)
        if uri.startswith(_base_uri):
            # Strip the base URI part
            uri = uri[len(_base_uri):]
        assert uri[0] == '#', "Annotation found on non-local URI"
        var_id = uri[1:] # Strip '#'
        var_objs = cellml_model.xml_xpath(u'*/cml:variable[@cmeta:id="%s"]' % var_id)
        assert len(var_objs) == 1, "Didn't find a single variable with ID " + var_id
        vars.append(var_objs[0])
    return vars
