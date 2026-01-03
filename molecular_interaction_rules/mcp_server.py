#!/usr/bin/env python3
#
# Molecular Interation Rules: MCP Server
# ---------------------------------------

# Imports
# -------

import json
import sys
from molecular_interaction_rules.molecular_database import MoleculerDatabase

class MCPServer(object):

    def __init__(self):
        
        self.database = MoleculerDatabase()
        self.tools = self._register_tools()

    def _register_tools(self):
        
        '''
        
        Register MCP Tools
        
        '''
        
        tools = [
            {
                'name': 'list_molecules',
                'description': 'Get a list of all available molecules in the database',
                'inputSchema': {
                    'type': 'object',
                    'properties': {},
                    'required': []
                }
            },
            {
                'name': 'get_molecule_fg_family',
                'description': 'Get the functional group family of a molecule',
                'inputSchema': {
                    'type': 'object',
                    'properties': {
                        'molecule': {
                            'type': 'string',
                            'description': 'Name of the molecule (e.g., "benzene")'
                        }
                    },
                    'required': ['molecule']
                }
            },
            {
                'name': 'get_atom_names',
                'description': 'Get the site names (atom names) in a molecule for interaction with another molecule',
                'inputSchema': {
                    'type': 'object',
                    'properties': {
                        'molecule': {
                            'type': 'string',
                            'description': 'Name of the molecule'
                        }
                    },
                    'required': ['molecule']
                }
            },
            {
                'name': 'get_monomer_coordinates',
                'description': 'Get the z-matrix coordinates for a monomer at a specific interaction site',
                'inputSchema': {
                    'type': 'object',
                    'properties': {
                        'molecule': {
                            'type': 'string',
                            'description': 'Name of the molecule'
                        },
                        'site': {
                            'type': 'string',
                            'description': 'Interaction site name (e.g., "RC1", "H1")'
                        }
                    },
                    'required': ['molecule', 'site']
                }
            },
            {
                'name': 'form_dimer_coordinates',
                'description': 'Generate z-matrix coordinates for a dimer formed from two monomers at specified interaction sites',
                'inputSchema': {
                    'type': 'object',
                    'properties': {
                        'molecule_1': {
                            'type': 'string',
                            'description': 'Name of first molecule'
                        },
                        'site_1': {
                            'type': 'string',
                            'description': 'Interaction site on first molecule'
                        },
                        'molecule_2': {
                            'type': 'string',
                            'description': 'Name of second molecule'
                        },
                        'site_2': {
                            'type': 'string',
                            'description': 'Interaction site on second molecule'
                        }
                    },
                    'required': ['molecule_1', 'site_1', 'molecule_2', 'site_2']
                }
            }
        ]
        
        return tools

    def handle_tool_call(self, tool_name, arguments):
        
        '''
        
        Handle Tool Execution
        
        '''
        
        try:
            if tool_name == 'list_molecules':
                result = self.database.get_molecule_list()
                return {'content': [{'type': 'text', 'text': json.dumps(result, indent=2)}]}
            
            elif tool_name == 'get_molecule_fg_family':
                molecule = arguments.get('molecule')
                result = self.database.get_molecule_fg_family(molecule)
                return {'content': [{'type': 'text', 'text': str(result)}]}
            
            elif tool_name == 'get_atom_names':
                molecule = arguments.get('molecule')
                result = self.database.get_atom_names(molecule)
                return {'content': [{'type': 'text', 'text': json.dumps(result, indent=2)}]}
            
            elif tool_name == 'get_monomer_coordinates':
                molecule = arguments.get('molecule')
                site = arguments.get('site')
                result = self.database.get_monomer_coordinates(molecule, site)
                return {'content': [{'type': 'text', 'text': str(result)}]}
            
            elif tool_name == 'form_dimer_coordinates':
                molecule_1 = arguments.get('molecule_1')
                site_1 = arguments.get('site_1')
                molecule_2 = arguments.get('molecule_2')
                site_2 = arguments.get('site_2')
                result = self.database.form_dimer_coordinates(
                    molecule_1, site_1, molecule_2, site_2
                )
                return {'content': [{'type': 'text', 'text': str(result)}]}
            
            else:
                return {'error': {'code': -32601, 'message': 'Tool not found'}}
        
        except Exception as e:
            return {'error': {'code': -32603, 'message': str(e)}}

    def handle_request(self, request):
        
        '''
        
        Handle MCP Request
        
        '''
        
        method = request.get('method')
        params = request.get('params', {})
        request_id = request.get('id')
        
        if method == 'tools/list':
            response = {
                'jsonrpc': '2.0',
                'id': request_id,
                'result': {'tools': self.tools}
            }
        
        elif method == 'tools/call':
            tool_name = params.get('name')
            arguments = params.get('arguments', {})
            result = self.handle_tool_call(tool_name, arguments)
            response = {
                'jsonrpc': '2.0',
                'id': request_id,
                'result': result
            }
        
        elif method == 'initialize':
            response = {
                'jsonrpc': '2.0',
                'id': request_id,
                'result': {
                    'protocolVersion': '2024-11-05',
                    'serverInfo': {
                        'name': 'molecular-interaction-rules-mcp',
                        'version': '1.0.0'
                    },
                    'capabilities': {
                        'tools': {}
                    }
                }
            }
        
        else:
            response = {
                'jsonrpc': '2.0',
                'id': request_id,
                'error': {'code': -32601, 'message': 'Method not found'}
            }
        
        return response

    def run(self):
        
        '''
        
        Run MCP Server via STDIO
        
        '''
        
        for line in sys.stdin:
            try:
                request = json.loads(line)
                response = self.handle_request(request)
                print(json.dumps(response), flush=True)
            except json.JSONDecodeError:
                error_response = {
                    'jsonrpc': '2.0',
                    'error': {'code': -32700, 'message': 'Parse error'}
                }
                print(json.dumps(error_response), flush=True)
            except Exception as e:
                error_response = {
                    'jsonrpc': '2.0',
                    'error': {'code': -32603, 'message': str(e)}
                }
                print(json.dumps(error_response), flush=True)

def main():
    
    '''
    
    MCP Server Entry Point
    
    '''
    
    server = MCPServer()
    server.run()

if __name__ == '__main__':
    main()
