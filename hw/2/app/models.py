from sqlalchemy import Column, String
from config import Base

class Molecule(Base):
    __tablename__ = "molecules"
    identifier = Column(String, primary_key=True, index=True)
    smile = Column(String)
