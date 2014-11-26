var utils = {

	uniqid: (function() {	
	    
		var i = 0;
	    
		return function(prefix) {       
	        var v = (new Date()).getTime();
	        return (prefix || "id") + v.toString() + i++;
	    }; 
	    
	})(),

	isLittleEndian: (function () {
		var buffer = new ArrayBuffer(2);
		new DataView(buffer).setInt16(0, 256, true);
		return new Int16Array(buffer)[0] === 256;
	})(),

	createElement: function (tag, attributes) {
		var el = document.createElement(tag);

		for (var name in attributes) {
			var att = document.createAttribute(name);
				att.value = attributes[name];

			el.setAttributeNode(att);
		}

		return el;

	}

};