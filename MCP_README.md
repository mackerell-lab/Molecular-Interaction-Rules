# Molecular Interaction Rules MCP Server

This package includes a Model Context Protocol (MCP) server that exposes the molecular interaction rules database for use with MCP-compatible AI assistants.

## Installation

Install the package with MCP support:

```bash
pip install molecular-interaction-rules
```

## Running the MCP Server

Run the MCP server using the command-line interface:

```bash
molecular-interaction-rules-mcp
```

The server runs via STDIO and communicates using JSON-RPC 2.0 protocol.

## MCP Configuration

To register this server with an MCP client, add the following to your MCP configuration file:

```json
{
  "mcpServers": {
    "molecular-interaction-rules": {
      "command": "molecular-interaction-rules-mcp"
    }
  }
}
```

## Available Tools

The MCP server provides the following tools:

### 1. list_molecules

Get a list of all available molecules in the database.

**Parameters:** None

**Example:**
```json
{
  "name": "list_molecules",
  "arguments": {}
}
```

### 2. get_molecule_fg_family

Get the functional group family of a molecule.

**Parameters:**
- `molecule` (string): Name of the molecule (e.g., "benzene")

**Example:**
```json
{
  "name": "get_molecule_fg_family",
  "arguments": {
    "molecule": "benzene"
  }
}
```

### 3. get_atom_names

Get the site names (atom names) in a molecule for interaction with another molecule.

**Parameters:**
- `molecule` (string): Name of the molecule

**Example:**
```json
{
  "name": "get_atom_names",
  "arguments": {
    "molecule": "benzene"
  }
}
```

### 4. get_monomer_coordinates

Get the z-matrix coordinates for a monomer at a specific interaction site.

**Parameters:**
- `molecule` (string): Name of the molecule
- `site` (string): Interaction site name (e.g., "RC1", "H1")

**Example:**
```json
{
  "name": "get_monomer_coordinates",
  "arguments": {
    "molecule": "benzene",
    "site": "RC1"
  }
}
```

### 5. form_dimer_coordinates

Generate z-matrix coordinates for a dimer formed from two monomers at specified interaction sites.

**Parameters:**
- `molecule_1` (string): Name of first molecule
- `site_1` (string): Interaction site on first molecule
- `molecule_2` (string): Name of second molecule
- `site_2` (string): Interaction site on second molecule

**Example:**
```json
{
  "name": "form_dimer_coordinates",
  "arguments": {
    "molecule_1": "benzene",
    "site_1": "RC1",
    "molecule_2": "benzene",
    "site_2": "RC1"
  }
}
```

## Direct Python Usage

You can also use the MCP server directly in Python:

```python
from molecular_interaction_rules import MCPServer

server = MCPServer()

# Handle a request
request = {
    'jsonrpc': '2.0',
    'id': 1,
    'method': 'tools/list',
    'params': {}
}

response = server.handle_request(request)
print(response)
```
