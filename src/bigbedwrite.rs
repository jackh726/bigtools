use crate::bbiwrite::{BBIWriteOptions};


pub struct BigBedWrite {
    pub path: String,
    pub options: BBIWriteOptions,
}

impl BigBedWrite {
    pub fn create_file(path: String) -> Self {
        BigBedWrite {
            path,
            options: BBIWriteOptions {
                compress: true,
                items_per_slot: 1024,
                block_size: 256,
            }
        }
    }
}