#!/usr/bin/env python3
#
# Molecular Interaction Rules: MCP CLI
# -------------------------------------

# Imports
# -------

import sys
import argparse
from molecular_interaction_rules.mcp_server import MCPServer

def main():
    
    '''
    
    CLI Entry Point for MCP Server
    
    '''
    
    parser = argparse.ArgumentParser(
        description='Molecular Interaction Rules MCP Server'
    )
    
    parser.add_argument(
        '--version',
        action='version',
        version='1.0.0'
    )
    
    args = parser.parse_args()
    
    server = MCPServer()
    server.run()

if __name__ == '__main__':
    main()
