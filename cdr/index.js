'use strict';

var path = require('path');

const CdrPlugin = require(path.join(global.CDR_INSTALL_DIR, '/commander/classes/CdrPlugin')).CdrPlugin;

module.exports = class CfeCdrPlugin extends CdrPlugin {
  constructor(urlBase) {
    super('pe', path.join(__dirname, 'web', urlBase));
  }

  getContent() {
    var result = {
      shortDescription: 'Position Estimator',
      longDescription: 'Position Estimator',
      nodes: {
				main: {
					type: CdrPlugin.ContentType.LAYOUT,
					shortDescription: 'Main',
					longDescription: 'Main PE.',
					filePath: '/main_layout.lyt'
				},
				hk: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE HK',
					longDescription: 'Position Estimator Housekeeping',
					filePath: '/hk.pug'
				},
				gps: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE GPS',
					longDescription: 'Position Estimator GPS',
					filePath: '/gps.pug'
				},
				baro: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE Baro',
					longDescription: 'Position Estimator Barometric Altimeter',
					filePath: '/baro.pug'
				},
				dist: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE Distance',
					longDescription: 'Position Estimator Distance Sensor',
					filePath: '/dist.pug'
				},
				land: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE Landingr',
					longDescription: 'Position Estimator Landing Sensor',
					filePath: '/land.pug'
				},
				flow: {
					type: CdrPlugin.ContentType.PANEL,
					shortDescription: 'PE Flow',
					longDescription: 'Position Estimator Flow Sensor',
					filePath: '/flow.pug'
				}
      }
    };

    return result;
  }
};

//var plugin = new CfeCdrPlugin();

//let CfeCdrPlugin = new class extends CdrPlugin {
