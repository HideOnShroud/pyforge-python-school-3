from typing import List, Optional, Generic, TypeVar
from pydantic import BaseModel, Field, ConfigDict


T = TypeVar('T')
class MoleculeSchema(BaseModel):
    identifier: str
    smile: str
    
    model_config = ConfigDict(from_attributes=True)



class RequestMolecule(BaseModel):
    parameters: MoleculeSchema = Field(...)
    
class Response(BaseModel, Generic[T]):
    code: str
    status: str
    message: str
    result: Optional[T]
    
